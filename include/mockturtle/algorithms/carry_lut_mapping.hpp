/*
  carry_carry_lut_mapping.hpp

  separate lists for carry nodes and carry lut nodes
  two method to check if node belongs to either carry or carry lut
*/

#pragma once

#include <cstdint>

#include <fmt/format.h>

#include "../utils/stopwatch.hpp"
#include "../views/topo_view.hpp"
#include "../networks/mig.hpp"
#include "cut_enumeration.hpp"
#include "cut_enumeration/mf_cut.hpp"
#include <kitty/print.hpp>

#define baseline 0
#define COST 4// 4// 4

namespace mockturtle
{

/*! \brief Parameters for carry_lut_mapping.
 *
 * The data structure `carry_lut_mapping_params` holds configurable parameters
 * with default arguments for `carry_lut_mapping`.
 */
struct carry_lut_mapping_params
{
  carry_lut_mapping_params()
  {
    cut_enumeration_ps.cut_size = 6;
    cut_enumeration_ps.cut_limit = 100;
  }

  /*! \brief Parameters for cut enumeration
   *
   * The default cut size is 6, the default cut limit is 8.
   */
  cut_enumeration_params cut_enumeration_ps{};

  /*! \brief Number of rounds for area flow optimization.
   *
   * The first round is used for delay optimization.
   */
  uint32_t rounds{3u};

  /*! \brief Number of rounds for exact area optimization. */
  uint32_t rounds_ela{1u};

  /*! \brief Be verbose. */
  bool verbose{true};

  /*! \brief Verbosity Level. >3 means print everything*/
  uint32_t verbosity = 6;

  /*! \brief Map to carry. */
  bool carry_mapping{true};

  /*! \brief Number of rounds for carry chain synthesis. */
  uint32_t max_rounds_carry{1000u};  
  
  /*! \brief Determines whether to print carry node combined LUT fcn. */
  bool carry_lut_combined{false};  
};

/*! \brief Statistics for carry_lut_mapping.
 *
 * The data structure `carry_lut_mapping_stats` provides data collected by running
 * `carry_lut_mapping`.
 */
struct carry_lut_mapping_stats
{
  /*! \brief Total runtime. */
  stopwatch<>::duration time_total{0};

  void report() const
  {
    std::cout << fmt::format( "[i] total time = {:>5.2f} secs\n", to_seconds( time_total ) );
  }
};

/* function to update all cuts after cut enumeration */
template<typename CutData>
struct carry_lut_mapping_update_cuts
{
  template<typename NetworkCuts, typename Ntk>
  static void apply( NetworkCuts const& cuts, Ntk const& ntk )
  {
    (void)cuts;
    (void)ntk;
  }
};

namespace detail
{

template<class Ntk, bool StoreFunction, typename CutData>
class carry_lut_mapping_impl
{
public:
  using network_cuts_t = network_cuts<Ntk, StoreFunction, CutData>;
  using cut_t = typename network_cuts_t::cut_t;

public:
  carry_lut_mapping_impl( Ntk& ntk, carry_lut_mapping_params const& ps, carry_lut_mapping_stats& st )
      : ntk( ntk ),
        ps( ps ),
        st( st ),
        flow_refs( ntk.size() ),
        map_refs( ntk.size(), 0 ),
        flows( ntk.size() ),
        delays( ntk.size() ),
        carry_cut_list( ntk.size() ),
        carry_nodes( ntk.size(), 0),
        carry_driver_nodes( ntk.size(), 0),
        carry_paths( ),
        cuts( cut_enumeration<Ntk, StoreFunction, CutData>( ntk, ps.cut_enumeration_ps ) )
  {
    carry_lut_mapping_update_cuts<CutData>().apply( cuts, ntk );
  }

  void run()
  {
    stopwatch t( st.time_total );

    /* compute and save topological order */
    top_order.reserve( ntk.size() );
    topo_view<Ntk>( ntk ).foreach_node( [this]( auto n ) {
      top_order.push_back( n );
    } );

    init_nodes();

    set_mapping_refs<false>();
    while ( iteration < ps.rounds )
    {
      compute_mapping<false>(false);
    }

    /* map to carry chain if true */
    if (ps.carry_mapping)
    {
      carry_chain_mapping();
    }

    while ( iteration < ps.rounds + ps.rounds_ela )
    {
      compute_mapping<true>(false);
    }
    derive_mapping();
  }

private:

  void carry_chain_mapping () {

      // Initialize carry chain mapping struct
      init_carry_chain_mapping();

      // Determine how many rounds of carry chain mapping
      uint32_t rounds_carry_mapping = 10;
      rounds_carry_mapping = num_worst_paths();
      if (rounds_carry_mapping > ps.max_rounds_carry)
        rounds_carry_mapping = ps.max_rounds_carry;
  
      // Iteratively map paths to carry chain
      for (uint32_t i = 0; i < rounds_carry_mapping; i++) { 

        if (ps.verbose) std::cout << "Carry Mapping Iteration " \
          << i << " targetting delay of " << delay << ".\n";
        if (ps.verbose) print_state();

        // Select nodes to be placed on carry 
        std::vector<node<Ntk>> path_for_carry_chain;
        if (!path_selection(path_for_carry_chain))
          break;

        std::cout << "size of path is " << path_for_carry_chain.size() << "\n";

        // Compute carry LUT mapping 
        compute_carry_mapping<false>(path_for_carry_chain);
        path_for_carry_chain.clear();
      }

      //remove_inverter_for_carry_mapping();
      //set_carry_mapping_refs<false>();

      // Delete me
      if (1) {
        std::cout << "There were " << carry_paths.size() << " paths placed on carry chain\n"; 
        //num_worst_paths();
      }
  }

  ///////////////////////////////////////////////////////////////
  // Selecting Path to be Placed on Carry Chain
  ///////////////////////////////////////////////////////////////

  // Depth first search 
  // Find path from output of LUT to target input LUT
  bool add_find_longest_path_in_mapping (std::vector<node<Ntk>>& path_for_carry_chain, 
      uint32_t index, uint32_t source_index, uint32_t dest_index) {

    auto const& n = ntk.index_to_node(source_index);
    
    if (source_index == dest_index) {
      return true;
    }

    if (is_a_carry_node(ntk.index_to_node(index))) {
      for ( auto leaf : carry_cut_list[index] )
        if (source_index == leaf) return false;
    } else {
      for ( auto leaf : cuts.cuts(index)[0] )
        if (source_index == leaf) return false;
    }
 
    bool found = false; 
    for (uint32_t i = 0; i < ntk.fanin_size(n); i++) {
      auto const leaf = ntk.get_children(n,i);
      auto const leaf_index = ntk.node_to_index(leaf);
      if (!is_a_carry_node(leaf) && !ntk.is_constant(leaf))
        found = add_find_longest_path_in_mapping (path_for_carry_chain, index, leaf_index, dest_index);
      if (found) {
        path_for_carry_chain.push_back(leaf);
        carry_nodes[leaf] += 1;
        if (ps.verbose && ps.verbosity > 2) std::cout << "adding " << leaf << "\n";
        return true;
      }
    }
    return found;
  }

  // Keep going deeper to its fanin until depth is found
  bool find_deepest_LUT (std::vector<node<Ntk>>& path_for_carry_chain, uint32_t index) { 
    
    if (ps.verbose && ps.verbosity > 2) std::cout << "Index " << index << "(" << delays[index] << ") ";

    if (delays[index] == 0) {
      return true;
    } 

    bool deepest = false;

    uint32_t max_leaf = index;
    uint32_t max_delay = 0;
    float min_fanout{std::numeric_limits<float>::max()};

    if (is_a_carry_node(ntk.index_to_node(index))) {
      std::cout << "carry\n";
      for ( auto leaf : carry_cut_list[index] ) {
        if (!is_a_carry_node(leaf)) { 
          if (ps.verbose && ps.verbosity > 2) std::cout << "\t" << leaf << ":" << delays[leaf] << "," << flow_refs[leaf] << "," << ntk.fanout_size(leaf) <<"\n"; 
          if (max_delay < delays[leaf] || (max_delay == delays[leaf] && (ntk.fanout_size(leaf) < min_fanout))) {
            max_leaf = leaf;
            max_delay = delays[leaf];
            min_fanout = flows[leaf];
          }
        }
      }
    } else {
      std::cout << "reg\n";
      for ( auto leaf : cuts.cuts( index )[0] ){
        if (!is_a_carry_node(leaf)) { 
          if (ps.verbose && ps.verbosity > 2) std::cout << "\t" << leaf << ":" << delays[leaf] << "," << flow_refs[leaf] << "," << ntk.fanout_size(leaf) <<"\n"; 
          if (max_delay < delays[leaf] || (max_delay == delays[leaf] && (ntk.fanout_size(leaf) < min_fanout))) {
            max_leaf = leaf;
            max_delay = delays[leaf];
            min_fanout = flows[leaf];
          }
        }
      }
    }

    if (max_leaf != index) {
      deepest = find_deepest_LUT(path_for_carry_chain, max_leaf);
      if (deepest)  {
        ntk.clear_visited();
        if(add_find_longest_path_in_mapping(path_for_carry_chain, index, index, max_leaf)) {
        } else { //cannot put this in carry
          return false;
        } 
      }
    } //else assert(0);

    return deepest;
  }
  
  bool path_selection (std::vector<node<Ntk>>& path_for_carry_chain) {

    for ( auto it = top_order.rbegin(); it != top_order.rend(); ++it ) {
      const auto node = *it;
      const auto index = ntk.node_to_index(node);

      // Node with worst delay
      if (delays[index] == delay && delay != 0 && !ntk.is_pi(node)) {
        std::cout << "\nFound target index " << index << "(" << delays[index] << ")";

        // Try to place a path starting from this node
        if (find_deepest_LUT(path_for_carry_chain,index)) {
          if (!is_a_carry_node(index)) {
            std::cout << "adding " << index << "\n";
            path_for_carry_chain.push_back(index);
            carry_nodes[index] += 1;
          }
          // Only place one path at a time
          break;
        } else {
          path_for_carry_chain.clear();
        }
      }
    }
    if (path_for_carry_chain.empty()) return false;

    // Keep the path information in carry_paths
    carry_paths.push_back(path_for_carry_chain);

    // Update driver info for figuring out which node drives which carry
    // in case there are 2 nodes being driving a carry node
    for (uint32_t j = 1; j < path_for_carry_chain.size(); j++) {
      carry_driver_nodes[path_for_carry_chain[j]] = path_for_carry_chain[j-1];
    }

    return true;
  }


  ///////////////////////////////////////////////////////////////
  // Remove inverter and check
  ///////////////////////////////////////////////////////////////

  // If the node is a carry, check its carry child for inversion
  void check_inverter (void) {

    for (auto carry_path: carry_paths) {
      if (carry_path.empty()) continue;
      for (uint32_t carry_i = 1; carry_i < carry_path.size()-1; carry_i++) {

        auto n = carry_path[carry_i];
        auto n_carry = carry_path[carry_i-1];

        if (ps.verbose && ps.verbosity > 3) std::cout << "node " << n << ": ";
        for (uint32_t i = 0; i < ntk.fanin_size(n); i++) {

          node<Ntk> child_node = ntk.get_children(n,i);  
          if (ps.verbose && ps.verbosity > 3) std::cout << child_node << "(" << ntk.is_complemented_children(n,i) << ")";
          if(child_node == n_carry) {
            if (ps.verbose && ps.verbosity > 3) std::cout << "*";
            if (ntk.is_complemented_children(n,i)) { 
              //if (carry_i == carry_path.size()-1) std::cout << "^";
              assert(0);
            }
          }
          if (ps.verbose && ps.verbosity > 3) std::cout << " ";
        } 
        if (ps.verbose && ps.verbosity > 3) std::cout << "\n";
      }
    } 
  }

  void remove_inverter_for_carry_mapping ( ) {
      // Must update cut truth tables since inverters were moved
      // Remove inverters in its path
      init_carry_chain_mapping();
      remove_inverter();
      cuts = cut_enumeration<Ntk, StoreFunction, CutData>( ntk, ps.cut_enumeration_ps ); 
      for (auto update_carry_path: carry_paths) {
        if (update_carry_path.empty()) break;
        compute_carry_mapping<true>(update_carry_path);
      }
      check_inverter();
  }

  void remove_inverter() {
    for (auto carry_path: carry_paths) {
      for (uint32_t carry_i = 1; carry_i < carry_path.size(); carry_i++) {

        auto carry_node = carry_path[carry_i];
        auto carry_child_node = carry_path[carry_i-1];
        
        bool complemented_carry_edge = false;
        for (uint32_t i = 0; i < ntk.fanin_size(carry_node); i++) {
          node<Ntk> child_node = ntk.get_children(carry_node,i);  
          if (child_node == carry_child_node && ntk.is_complemented_children(carry_node,i))
            complemented_carry_edge = true;
        }

        // If the carry path is complemented, go through all the nodes
        // and flip its children that match this node and the node itself
        if (complemented_carry_edge) {

          // flip children
          for (uint32_t i = 0; i < ntk.fanin_size(carry_node); i++) {
            node<Ntk> child_node = ntk.get_children(carry_node,i);  
            if (ps.verbose && ps.verbosity > 3) std::cout << "\t\tflipping child " << child_node << "\n";
            ntk.flip_children(carry_node, i);
          }
          
          ntk.foreach_node( [&]( auto n, auto ) {
            for (uint32_t i = 0; i < ntk.fanin_size(n); i++) {
              if (ntk.get_children(n,i) == carry_node) {
                if (ps.verbose && ps.verbosity > 3) std::cout << "\t\tflipping child " << ntk.get_children(n,i) << \
                    " of node " << n << "\n";
                if (ps.verbose && ps.verbosity > 3) std::cout << "\t\t\tfrom " << ntk.is_complemented_children(n,i); 
                ntk.flip_children(n, i);
                if (ps.verbose && ps.verbosity > 3) std::cout << " to " << ntk.is_complemented_children(n,i) << "\n"; 
              }
            }
          });
          ntk.foreach_po( [&]( auto const& s ) {
            if (ntk.get_node(s) == carry_node) {
              if (ps.verbose && ps.verbosity > 3)std::cout << "\t\tflipping output " << ntk.get_node(s) << "\n";
              ntk.flip_complement_output(s);
            }
          });

        }
      }
    }
  }

  ///////////////////////////////////////////////////////////////
  // Initialization
  ///////////////////////////////////////////////////////////////

  void init_carry_chain_mapping() {
    //do nothing for now
  }

  void init_nodes() {
    ntk.foreach_node( [this]( auto n, auto ) {
      const auto index = ntk.node_to_index( n );
      map_refs[index] = 0;

      if ( ntk.is_constant( n ) || ntk.is_pi( n ) ) {
        /* all terminals have flow 1.0 */
        flow_refs[index] = 1.0f;
      } else {
        flow_refs[index] = static_cast<float>( ntk.fanout_size( n ) );
      }

      flows[index] = cuts.cuts( index )[0]->data.flow;
      delays[index] = cuts.cuts( index )[0]->data.delay*LUT_DELAY;
    } );
  }

  ///////////////////////////////////////////////////////////////
  // Map to LUTs before carry nodes 
  ///////////////////////////////////////////////////////////////

  template<bool SetLUT>
  void compute_carry_mapping (std::vector<node<Ntk>> path_for_carry_chain)
  {

    node<Ntk> n_carryin, n_first, n_second;

    // This function should not be called with an empty carry path
    assert(!path_for_carry_chain.empty());

    // Map nodes in pairs (i+=2) since each ALM has 2 adders
    for (uint32_t i = 1; i < path_for_carry_chain.size(); i+=2) {

      // Assign carryin, 1st/2nd node
      n_carryin = path_for_carry_chain[i-1];
      n_first = path_for_carry_chain[i];
      n_second = (i == path_for_carry_chain.size()-1) ? 0 : path_for_carry_chain[i+1];

      if (ps.verbose && ps.verbosity > 3) {
        std::cout << "i " << i << ": " << n_carryin << " " << n_first << " " << n_second << "\n";
      }

      compute_carry_LUT_mapping<SetLUT>(n_first, n_second, n_carryin);
      //select_no_cuts<SetLUT>(n_first, n_second, n_carryin);
    }
  }

  // This function takes 2 index with its carry-in index and figures out 
  // if it cuts can be mapped to LUTs before carry node and the cost
  template <bool SetLUT>
  int compute_carry_LUT_mapping(uint32_t i1, uint32_t i2, uint32_t ic) {

    double cost = 0;
    double max_cost = -1;
    int32_t max_i1_1 = -1;
    int32_t max_i1_2 = -1;
    int32_t max_i2_1 = -1;
    int32_t max_i2_2 = -1;

    // Get relevant children node for mapping
    node<Ntk> i1_child[2] = {0};
    node<Ntk> i2_child[2] = {0};
    get_children_node (i1_child, i1, ic); 
    if (i2 != 0) get_children_node (i2_child, i2, i1); 

    uint32_t total_num_cuts = 0;
    uint32_t total_num_legal_cuts = 0;

    for (uint32_t cut_i1_1 = 0; cut_i1_1 < cuts.cuts(i1_child[0]).size(); cut_i1_1++) {
      for (uint32_t cut_i1_2 = 0; cut_i1_2 < cuts.cuts(i1_child[1]).size(); cut_i1_2++) {
        for (uint32_t cut_i2_1 = 0; cut_i2_1 < cuts.cuts(i2_child[0]).size(); cut_i2_1++) {
          for (uint32_t cut_i2_2 = 0; cut_i2_2 < cuts.cuts(i2_child[1]).size(); cut_i2_2++) {

        total_num_cuts++;

        // Array holding the nodes to the halfs of ALM as separate 4-LUT
        // [0][0] and [1][0] contains # of input
        uint32_t i1_child_cuts[2][5] = {0}; 
        uint32_t i2_child_cuts[2][5] = {0}; 

        uint32_t abce0f0[5] = {0};
        uint32_t abde1f1[5] = {0};
  
        uint32_t ntotal_i1 = 0; uint32_t ntotal_i2 = 0;
        uint32_t nshared_i1 = 0; uint32_t nshared_i2 = 0;

        bool newcheck = check_cut_legality(i1_child[0], i1_child[1], i2_child[0], i2_child[1], \
          cut_i1_1, cut_i1_2, cut_i2_1, cut_i2_2);

        bool oldcheck = true;          
        ////////////////////////////////
        // Checking input requirements
        //////////////////////////////// 
        if (insert_cuts_to_array (i1_child_cuts, i1, ic, cut_i1_1, cut_i1_2)) {
          if (!check_5lut_legality (i1_child_cuts, abce0f0, &nshared_i1, &ntotal_i1))
            oldcheck = false;
        } else oldcheck = false;

        // if i2 is 0, it's not used 
        if (i2 != 0) {
          if (insert_cuts_to_array (i2_child_cuts, i2, i1, cut_i2_1, cut_i2_2)) { 
            if (!check_5lut_legality (i2_child_cuts, abde1f1, &nshared_i2, &ntotal_i2))
              oldcheck = false;
          } else oldcheck = false;
        } 

        uint32_t nshared = 0; uint32_t ntotal = 0;
        if (!check_6lut_legality(abce0f0, abde1f1, nshared_i1, nshared_i2, \
          ntotal_i1, ntotal_i2, nshared, ntotal)) {
          oldcheck = false;
        }

        // TODO: when new check is finished
        if (!oldcheck) continue;
        total_num_legal_cuts++;

        ////////////////////////////////
        // Calculating cost of cut 
        //////////////////////////////// 
        double cost = cost_of_cuts<COST>(i1, i2, i1_child[0], i1_child[1],
            i2_child[0], i2_child[1], cut_i1_1, cut_i1_2, 
            cut_i2_1, cut_i2_2, nshared_i1, nshared_i2, 
            ntotal_i1, ntotal_i2, nshared);

        if (cost > max_cost) {
          max_cost = cost;
          max_i1_1 = cut_i1_1;
          max_i1_2 = cut_i1_2;
          max_i2_1 = cut_i2_1;
          max_i2_2 = cut_i2_2;
        }
        if (ps.verbose && ps.verbosity > 3) {
          std::cout << "Cut " << cut_i1_1 << " " << cut_i1_2 << " " << cut_i2_1 << " " << cut_i2_2 << ": ";
          std::cout << cost << "\n";
          std::cout << "\t" << cuts.cuts(i1_child[0])[cut_i1_1] << " " << cuts.cuts(i1_child[0])[cut_i1_1]->data.delay*20 << "\n";
          std::cout << "\t" << cuts.cuts(i1_child[1])[cut_i1_2] << " " << cuts.cuts(i1_child[1])[cut_i1_2]->data.delay*20 << "\n";
          std::cout << "\t" << cuts.cuts(i2_child[0])[cut_i2_1] << " " << cuts.cuts(i2_child[0])[cut_i2_1]->data.delay*20 << "\n";
          std::cout << "\t" << cuts.cuts(i2_child[1])[cut_i2_2] << " " << cuts.cuts(i2_child[1])[cut_i2_2]->data.delay*20 << "\n";
        }
          }
        }
      }
    }

    if (ps.verbose && ps.verbosity > 3) {
      std::cout << "Total number of cuts are " << total_num_cuts << " and " << total_num_legal_cuts << " were legal\n";
      std::cout << "Selected cut is " << max_i1_1 << " " << max_i1_2 << " " << max_i2_1 << " " <<  max_i2_2 << "\n";
    }

    insert_cut_to_carry_list<SetLUT>(i1, ic, i1_child[0], i1_child[1], max_i1_1, max_i1_2); 
    insert_cut_to_carry_list<SetLUT>(i2, i1, i2_child[0], i2_child[1], max_i2_1, max_i2_2); 

    //update_mapping_ref();
/*     
    node<Ntk> slowest_leaf;
    uint32_t time = 0; 
    if (ps.verbose && ps.verbosity > 3) std::cout << "Node* " << ntk.node_to_index(index1) << ":";
      for ( auto leaf : carry_cut_list[index1] ) {
        if (delays[leaf] > time) {
          time = delays[leaf];
          slowest_leaf = leaf;
        }
        if (ps.verbose && ps.verbosity > 3)std::cout << leaf << "(" << delays[leaf] << ") ";
    }
    std::cout << "\n";
    if (carry_driver_nodes[index1] == slowest_leaf)
      time += CARRY_DELAY;
    else 
      time += LUT_ADDER_DELAY;
    delays[index1] = time;

    time = 0;
    if (ps.verbose && ps.verbosity > 3) std::cout << "Node* " << ntk.node_to_index(index2) << ":";
      for ( auto leaf : carry_cut_list[index2] ) {
        if (delays[leaf] > time) {
          time = delays[leaf];
          slowest_leaf = leaf;
        }
        if (ps.verbose && ps.verbosity > 3)std::cout << leaf << "(" << delays[leaf] << ") ";
    }
    std::cout << "\n";
    if (carry_driver_nodes[index2] == slowest_leaf)
      time += CARRY_DELAY;
    else 
      time += LUT_ADDER_DELAY;
    delays[index2] = time;
*/
    return 0;
  } 

  // This function takes in an input array of both LUTs
  // It looks for the unique set of inputs
  uint32_t find_unique_set (uint32_t* n_shared, uint32_t input_array[2][5], \
      uint32_t unique_array[8]) {

    bool match = false;
    uint32_t index = 0;

    for (uint32_t i = 0; i < input_array[0][0]; i++) {
      for (uint32_t j = 0; j < input_array[1][0]; j++) {
        if (input_array[0][i+1] == input_array[1][j+1]) {
          unique_array[*n_shared] = input_array[0][i+1];  
          (*n_shared)++;
          index++;
          match = true;
          continue;
        } 
      }
      if (!match) {
        unique_array[*n_shared] = input_array[0][i+1];  
          (*n_shared)++;
        index++;
      }
      match = false;
    }
    for (uint32_t i = 0; i < input_array[1][0]; i++) {
      for (uint32_t j = 0; j < input_array[0][0]; j++) {
        if (input_array[0][i+1] == input_array[1][j+1]) {
          match = true;
          continue;
        } 
      }
      if (!match) {
        unique_array[index] = input_array[1][i+1];  
          (*n_shared)++;
        index++;
      }
      match = false;
    } 

    return index;
  }

  bool check_5lut_legality (uint32_t input_array[2][5], \
      uint32_t unique_array[5], uint32_t* n_shared, uint32_t* n_total) {

    bool match = false;

    uint32_t abc = 0;
    uint32_t e0f0 = 0;

    // Find all shared values and put it in first 3 of the unique_array
    // The unique values from array[0] (doesn't match) go into the last 2 elements of array
    for (uint32_t i = 0; i < input_array[0][0]; i++) {
      for (uint32_t j = 0; j < input_array[1][0]; j++) {
        if (input_array[0][i+1] == input_array[1][j+1]) {
          unique_array[abc] = input_array[0][i+1];  
          abc++;
          match = true;
          continue;
        } 
      }
      if (!match) {
        e0f0++;
        if (e0f0 > 5) return false;
        unique_array[5-e0f0] = input_array[0][i+1];  
      }
      match = false;
    }
    // The unique values from arrayp1[ (doesn't match) go into the last 2 elements of array
    match = false;
    for (uint32_t i = 0; i < input_array[1][0]; i++) {
      for (uint32_t j = 0; j < input_array[0][0]; j++) {
        if (input_array[0][j+1] == input_array[1][i+1]) {
          match = true;
          continue;
        } 
      }
      if (!match) {
        e0f0++;
        if (e0f0 > 5) return false;
        unique_array[5-e0f0] = input_array[1][i+1];  
      }
      match = false;
    } 

    // abc + e0f0 should be less than 5
    // n_shared keeps track of how many are actually shared 
    (*n_total) = abc + e0f0;
    (*n_shared) = abc;
    if (abc <= 3 && abc + e0f0 <= 5) return true;
    return false;
  }

  bool check_6lut_legality (uint32_t abce0f0[5], uint32_t abde1f1[5], \
      uint32_t nshared_index1, uint32_t nshared_index2, \
      uint32_t ntotal_index1, uint32_t ntotal_index2, 
      uint32_t& nshared_matches, uint32_t& nmatches) {

    // At this point, both 5 LUTs are legal
    // Now we need to check for 6 LUT input legality
    // Depending on the total number of inputs, 

    // From the shared list, how many shared inputs (for a or b)
    for (uint32_t i = 0; i < nshared_index1; i++) {
      assert (abce0f0[i] != 0);
      for (uint32_t j = 0; j < nshared_index2; j++) {
        if (abce0f0[i] == abde1f1[j]) nshared_matches++;
      }
    }
    // From all the list, how many shared inputs
    for (uint32_t i = 0; i < 5; i++) {
      if (abce0f0[i] == 0) continue;
      for (uint32_t j = 0; j < 5; j++) {
        if (abce0f0[i] == abde1f1[j]) nmatches++;
      }
    }
    // There are 8 pins available
    if ( (ntotal_index1 + ntotal_index2 - nmatches) > 8 )
      return false;

    // There only limited shared pins available
    if ( nshared_index1 > 3 || nshared_index2 > 3 )
      return false; 

 
    if ( ( ntotal_index1 > 3 ) && ((ntotal_index2 - nmatches) > 2) || \
         ( ntotal_index2 > 3 ) && ((ntotal_index1 - nmatches) > 2))
       return false;
    
    if (int(nshared_matches) >= int((nshared_index1 - 1) + (nshared_index2 - 1) - \
          (5 - ntotal_index1) - (5 - ntotal_index2))) { 
      return true;
    } 

    //assert(0);
    std::cout << "might be rejecting legal cut here\n";
    return false;
  }

  // This function adds child's inputs to an array for input legality checking
  // input_array [LUT number][0] holds # of values to follow
  void insert_child_to_array (uint32_t input_array[2][5], node<Ntk>child_index[2], 
      uint32_t index, uint32_t carryin) {

    uint32_t curr_LUT = 0;
    uint32_t curr_input= 0;
    auto const& curr_node = ntk.index_to_node(index);
    assert (ntk.fanin_size(curr_node) == 3);
    for (uint32_t i_fanin = 0; i_fanin < ntk.fanin_size(curr_node); i_fanin++) { 
      node<Ntk> leaf_node = ntk.get_children (curr_node, i_fanin);  
      
     
      if ((leaf_node == carryin) | (leaf_node == index)) {
        continue;
      } else {
        assert (curr_LUT < 2);

        if (ntk.is_constant(leaf_node)) {
          
        } else if (ntk.is_pi(leaf_node)) {
          input_array[curr_LUT][0]++;
          input_array[curr_LUT][curr_input+1] = leaf_node;
          curr_input++;
          child_index[curr_LUT] = leaf_node;
        } else {
          // Check the children of the cut for inputs
          for (uint32_t i_fanin = 0; i_fanin < ntk.fanin_size(leaf_node); i_fanin++) { 
            node<Ntk> child_leaf_node = ntk.get_children (leaf_node, i_fanin);  
            if (ntk.is_constant(child_leaf_node)) continue;
            input_array[curr_LUT][0]++;
            input_array[curr_LUT][curr_input+1] = child_leaf_node;
            curr_input++;
          }
          child_index[curr_LUT] = leaf_node;
        }
        curr_LUT++;
        curr_input = 0;
      }
    }
  }

  void get_children_node (node<Ntk>child_index[2], uint32_t index, uint32_t carryin) {

    uint32_t curr_LUT = 0;
    auto const& curr_node = ntk.index_to_node(index);

    // Why? 
    std::cout << curr_node << "(" << ntk.fanin_size(curr_node) << ")\n";
    assert (ntk.fanin_size(curr_node) == 3);

    for (uint32_t i_fanin = 0; i_fanin < ntk.fanin_size(curr_node); i_fanin++) { 
      node<Ntk> leaf_node = ntk.get_children (curr_node, i_fanin);  
      //std::cout << "\t" << leaf_node << "\n";
      
      // If the chilren node is the carry-in, skip
      if (leaf_node == carryin) {
        continue;
      } else {
        child_index[curr_LUT] = leaf_node;
        curr_LUT++;
      }
    }
  }

  // This function takes an array and an index with its carry-in index
  // and figures out the 2 driving children and fills input array with chilldren's cut nodes
  bool insert_cuts_to_array (uint32_t cut_array[2][5], 
      uint32_t index, uint32_t carryin, uint32_t cut_index1, uint32_t cut_index2) {

    uint32_t curr_LUT = 0;
    uint32_t curr_input= 0;
    auto const& curr_node = ntk.index_to_node(index);

    // This makes sure that we are dealing with MIG nodes
    assert (ntk.fanin_size(curr_node) == 3);

    for (uint32_t i_fanin = 0; i_fanin < ntk.fanin_size(curr_node); i_fanin++) { 
      node<Ntk> leaf_node = ntk.get_children (curr_node, i_fanin);  

      
      // If the chilren node is the carry-in, skip
      if (leaf_node == carryin) {
        continue;
      } else {
         
        // If the child node is PI, no cuts so add to it 
        if (ntk.is_pi(leaf_node)) {
          cut_array[curr_LUT][0]++;
          cut_array[curr_LUT][curr_input+1] = leaf_node;
          curr_input++;

        // If the child node is not a constant, check its cuts
        } else if (!ntk.is_constant(leaf_node)){

          // Check each cut of leaf node and add to the array
          auto cut_index = cut_index1;
          if (curr_LUT == 1) cut_index = cut_index2;
          cut_t const& cut = cuts.cuts(leaf_node)[cut_index]; 

          if (cut.size() > 4 || cut.size() == 0) {
            return false;
          }
          for (auto cut_leaf: cut) { 
            if (ntk.is_constant(cut_leaf)) continue;
            cut_array[curr_LUT][0]++;
            cut_array[curr_LUT][curr_input+1] = cut_leaf;
            curr_input++;
          }
        }
        curr_LUT++;
        curr_input = 0;
      }
    }
    return true;
  }
  uint32_t get_bit (uint32_t val, uint32_t pos) {
    return (val >> pos)&0x1;
  }
  uint32_t get_not_bit (uint32_t val, uint32_t pos) {
    if (((val >> pos)&0x1) == 1) return 0;
    return 0x1;
  }

  kitty::dynamic_truth_table compute_carry_function (bool child_complement[3], uint32_t unique_cuts[5], uint32_t n_total, 
      cut_t const& cut1, cut_t const& cut2) {

    // There cannot be more than 5 inputs to the 2 4-LUTs with the extra carry chain
    assert ( n_total + 1 <= 6);

    // First find the index order from the unique list
    uint32_t cut1_index[5] = {0};
    uint32_t cut2_index[5] = {0};

    uint32_t k = 0;
    for (auto cut_leaf: cut1) { 
      for (uint32_t i = 0; i < n_total; i++) {
        if (cut_leaf == unique_cuts[i]) {
          cut1_index[k] = i;
          k++;
        }
      }
    }
    k = 0;
    for (auto cut_leaf: cut2) { 
      for (uint32_t i = 0; i < n_total; i++) {
        if (cut_leaf == unique_cuts[i]) {
          cut2_index[k] = i;
          k++;
        }
      }
    }
    
    // Create a truth table
    uint32_t lut_input_size = n_total + 1; 
    kitty::dynamic_truth_table function (lut_input_size);
    function._bits[0]= 0x0000000000000000;

    auto cut1_function = cuts.truth_table( cut1 ); 
    auto cut2_function = cuts.truth_table( cut2 ); 
    
    // Go through each combination of the input
    // i4 i3 i2 i1 i0
    // 0 0 0 0 0
    // 0 0 0 0 1
    //std::cout << "\tFunction combination\n";
    uint32_t n_combinations = pow(2,n_total+1);
 
    // For actual functionality (CEC)
    if (ps.carry_lut_combined) {
      for (uint32_t i = 0; i < n_combinations; i++) {
        //  Figure out which bit from the LUT represents that combination
        uint32_t cut1_i = 0;
        for (uint32_t j = 0; j < cut1.size(); j++) {
          cut1_i += (get_bit(i, cut1_index[j]) << j);
        } 
        uint32_t cut2_i = 0;
        for (uint32_t j = 0; j < cut2.size(); j++) {
          cut2_i += (get_bit(i, cut2_index[j]) << j);
        } 

        //  With the bit, find the output of LUT
        uint32_t cut1_val = child_complement[0]? get_not_bit(cut1_function._bits[0], cut1_i) :
            get_bit(cut1_function._bits[0], cut1_i); 
        uint32_t cut2_val = child_complement[1]? get_not_bit(cut2_function._bits[0], cut2_i) :
            get_bit(cut2_function._bits[0], cut2_i); 
        uint32_t carry_val = child_complement[2]? get_not_bit (i, n_total) : get_bit(i, n_total);

        //  Combine those with the extra carry in with majority fcn
        uint64_t lut_val = ((cut1_val && cut2_val) || (cut1_val && carry_val) || (cut2_val & carry_val));
        
        //  Set that value in the right place of the function
        function._bits[0] += uint64_t( lut_val << i);
      }
    } else {
      for (uint32_t i = 0; i < n_combinations; i++) {
         //  Figure out which bit from the LUT represents that combination
        uint32_t cut1_i = 0;
        for (uint32_t j = 0; j < cut1.size(); j++) {
          cut1_i += (get_bit(i, cut1_index[j]) << j);
        } 
        uint32_t cut2_i = 0;
        for (uint32_t j = 0; j < cut2.size(); j++) {
          cut2_i += (get_bit(i, cut2_index[j]) << j);
        } 

        //  With the bit, find the output of LUT
        uint32_t cut1_val = child_complement[0]? get_not_bit(cut1_function._bits[0], cut1_i) :
            get_bit(cut1_function._bits[0], cut1_i); 
        uint32_t cut2_val = child_complement[1]? get_not_bit(cut2_function._bits[0], cut2_i) :
            get_bit(cut2_function._bits[0], cut2_i); 
        uint32_t carry_val = child_complement[2]? get_not_bit (i, n_total) : get_bit(i, n_total);

        //  Combine those with the extra carry in with majority fcn
        uint64_t lut_val = carry_val? cut1_val : cut2_val; 
        
        //  Set that value in the right place of the function
        function._bits[0] += uint64_t( lut_val << i);
      }
    }
    //std::cout << "\tFunction:"; 
    //kitty::print_hex(function);
    //std::cout << "\n";
    return function;
  }
  std::array<signal<Ntk>, 3> ordered_children( node<Ntk> const& n ) const
  {
    std::array<signal<Ntk>, 3> children;
    ntk.foreach_fanin( n, [&children]( auto const& f, auto i ) { children[i] = f; } );
    //std::sort( children.begin(), children.end(), [this]( auto const& c1, auto const& c2 ) {
    //  return ntk.level( ntk.get_node( c1 ) ) < ntk.level( ntk.get_node( c2 ) );
    //} );
    return children;
  }

  void remove_zero_array (uint32_t unique_array[5]) {

    uint32_t nonzero_i = 0;

    for (uint32_t i = 0; i < 5; i++) {
      if (unique_array[i] != 0) {
        unique_array[nonzero_i] = unique_array[i];
        if (i > nonzero_i)
          unique_array[i] = 0;
        nonzero_i++;
      }
    }
  }

  // convention is that carry_cut_list contains 
  // the inputs to the carry LUT (so push_back)*
  //  *maybe change to keep cut numbers
  template<bool SetLUT>
  void insert_cut_to_carry_list(uint32_t index, uint32_t cindex_c,
      uint32_t cindex_1, uint32_t cindex_2, int32_t i, int32_t j) {

    if (index == 0) return;

    auto n = ntk.index_to_node(index);

    bool child_complement[3] = {0};
    //std::cout << "children are ";
    ntk.foreach_fanin( n, [&]( auto const& f ) {
      if (cindex_1 == ntk.node_to_index(ntk.get_node(f))) {
        child_complement[0] = ntk.is_complemented(f);
      } else if (cindex_2 == ntk.node_to_index(ntk.get_node(f))) {
        child_complement[1] = ntk.is_complemented(f);
      } else if (cindex_c == ntk.node_to_index(ntk.get_node(f))) {
        child_complement[2] = ntk.is_complemented(f);
      } else { 
        assert(0 && "Children index do not match with fanin.");
      }
      //std::cout << ntk.get_node(f) << "(" << ntk.is_complemented(f) << "),";
    } );
    //std::cout << "\n";

    carry_cut_list[index].clear();

    if (i < 0 || j < 0) {
      carry_cut_list[index].push_back(cindex_1);
      carry_cut_list[index].push_back(cindex_2);
      carry_cut_list[index].push_back(cindex_c);

      if (SetLUT) {
        kitty::dynamic_truth_table mig_function (3);
        mig_function._bits[0]= 0x0000000000000000;
        for (uint32_t i = 0; i < 8; i++) {    
          uint32_t a = child_complement[0]? get_not_bit(i,0):get_bit(i,0); 
          uint32_t b = child_complement[1]? get_not_bit(i,1):get_bit(i,1); 
          uint32_t c = child_complement[2]? get_not_bit(i,2):get_bit(i,2);

          //  Combine those with the extra carry in with majority fcn
          uint64_t lut_val = ((a && b) || (b && c) || (a & c));
        
          mig_function._bits[0] += uint64_t( lut_val << i);
        }
        ntk.set_cell_function(n, mig_function );
      }
      return;
    }

    assert(i>=0 && j>=0);
    cut_t const& cut_1 = cuts.cuts(cindex_1)[i];
    cut_t const& cut_2 = cuts.cuts(cindex_2)[j];

    // Set LUT function
    // cut_1 and cut_2 with MIG carry node should be
    // turned into truth table 
    uint32_t index_child_cuts[2][5] = {0}; 
    uint32_t unique_cuts[5] = {0};
    uint32_t n_shared = 0;   
    uint32_t n_total = 0;   
    insert_cuts_to_array(index_child_cuts, index, cindex_c, i, j);
    check_5lut_legality(index_child_cuts, unique_cuts, &n_shared, &n_total);
    remove_zero_array(unique_cuts);


    if (SetLUT) {
      kitty::dynamic_truth_table function = compute_carry_function(child_complement, unique_cuts, n_total, cut_1, cut_2);
      ntk.set_cell_function(n, function);
    }

    for (uint32_t i = 0; i < n_total; i++) {
      carry_cut_list[index].push_back(unique_cuts[i]);
    }
    carry_cut_list[index].push_back(cindex_c);
  }
/*
  template <bool SetLUT>
  void set_function (node<Ntk> n, uint32_t cindex_1, uint32_t cindex_2, uint32_t cindex_c) {
    bool child_complement[3] = {0};
    //std::cout << "node " << n << " children are ";
    ntk.foreach_fanin( n, [&]( auto const& f ) {
      if (cindex_1 == ntk.node_to_index(ntk.get_node(f))) {
        child_complement[0] = ntk.is_complemented(f);
      } else if (cindex_2 == ntk.node_to_index(ntk.get_node(f))) {
        child_complement[1] = ntk.is_complemented(f);
      } else if (cindex_c == ntk.node_to_index(ntk.get_node(f))) {
        child_complement[2] = ntk.is_complemented(f);
      } else { 
        std::cout << "child is " << ntk.get_node(f) << " and does not fit " << \
            cindex_1 << " " << cindex_2 << " " << cindex_c << "\n";
        assert(0 && "Children index do not match with fanin.");
      }
    } );


    if (SetLUT) {
      kitty::dynamic_truth_table mig_function (3);
      mig_function._bits[0]= 0x0000000000000000;
      for (uint32_t i = 0; i < 8; i++) {    
        uint32_t a = child_complement[0]? get_not_bit(i,0):get_bit(i,0); 
        uint32_t b = child_complement[1]? get_not_bit(i,1):get_bit(i,1); 
        uint32_t c = child_complement[2]? get_not_bit(i,2):get_bit(i,2);

        //  Combine those with the extra carry in with majority fcn
        uint64_t lut_val = ((a && b) || (b && c) || (a & c));
      
        mig_function._bits[0] += uint64_t( lut_val << i);
      }
      ntk.set_cell_function(n, mig_function );
    }
 
  }

  template <bool SetLUT>
  void select_no_cuts(uint32_t index1, uint32_t index2, uint32_t indexc) {

    node<Ntk> index1_child[2] = {0};
    node<Ntk> index2_child[2] = {0};

    get_children_node (index1_child, index1, indexc); 
    carry_cut_list[index1].clear();
    carry_cut_list[index1].push_back(index1_child[0]); 
    carry_cut_list[index1].push_back(index1_child[1]); 
    carry_cut_list[index1].push_back(indexc); 
  
    set_function<SetLUT>(ntk.index_to_node(index1), index1_child[0], index1_child[1], indexc);

    if (index2 != 0) {
      get_children_node (index2_child, index2, index1); 
      carry_cut_list[index2].clear();
      carry_cut_list[index2].push_back(index2_child[0]); 
      carry_cut_list[index2].push_back(index2_child[1]); 
      carry_cut_list[index2].push_back(index1); 
      set_function<SetLUT>(ntk.index_to_node(index2), index2_child[0], index2_child[1], index1);
    }

  }
*/
  // Find the delay of the specific cut
  // It should look for the slowest of all nodes
  uint32_t delay_of_cut (uint32_t index, uint32_t cut_index) {
   uint32_t cut_delay = 0; 
      for (uint32_t cut: cuts.cuts(index)[cut_index]) {
        if (delays[cut] > cut_delay) {
          cut_delay = delays[cut];
        }
      }
    if (cut_delay == 0) return 1;
    else return cut_delay;
  }

  template <uint32_t cost_type> double cost_of_cuts(uint32_t index1, uint32_t index2,
      uint32_t index1_1, uint32_t index1_2,
      uint32_t index2_1, uint32_t index2_2, uint32_t i1, uint32_t i2, uint32_t j1, uint32_t j2,
      uint32_t nshared1, uint32_t nshared2, uint32_t ntotal1, uint32_t ntotal2, uint32_t nshared) {

    double cost = 0;
      
    if (cost_type == 0) {
      cost += cuts.cuts(index1_1)[i1].size() + cuts.cuts(index1_2)[i2].size();
      cost += cuts.cuts(index2_1)[j1].size() + cuts.cuts(index2_2)[j2].size();

    } else if (cost_type == 1) { // Best Delay
  
      cost += (1.0/(double)delay_of_cut (index1_1, i1));
      cost += (1.0/(double)delay_of_cut (index1_2, i2));
      cost += (1.0/(double)delay_of_cut (index2_1, j1));
      cost += (1.0/(double)delay_of_cut (index2_2, j2));

    } else if (cost_type == 2) { // Most Number of Shared Inputs
  
      cost += 2*nshared; 
      cost += nshared1;
      cost += nshared2;
    
    } else if (cost_type == 3) { // Delay with shared input

      cost += (1/(double)delay_of_cut (index1_1, i1));
      cost += (1/(double)delay_of_cut (index1_2, i2));
      cost += (1/(double)delay_of_cut (index2_1, j1));
      cost += (1/(double)delay_of_cut (index2_2, j2));

      cost += 0.5*(nshared + nshared + nshared2); 

    } else if (cost_type == 4) {
      cost += (1/(double)delay_of_cut (index1_1, i1));
      cost += (1/(double)delay_of_cut (index1_2, i2));
      cost += (1/(double)delay_of_cut (index2_1, j1));
      cost += (1/(double)delay_of_cut (index2_2, j2));
      cost *= 20; 
      cost *= (1 + nshared1 + nshared2 + nshared);
      
    } else {
      assert(0);
    }
    return cost;
  }

  uint32_t number_of_match_in_cuts (auto cut_1, auto cut_2) {
    uint32_t match = 0;
    for (auto leaf_1: cut_1) {
      for (auto leaf_2: cut_2) {
        if (leaf_1 == leaf_2) match++;
      }
    }
    return match;
  }

  void match_vector (std::vector<node<Ntk>> return_vector, auto cut_1, auto cut_2) {
    uint32_t match = 0;
    for (auto leaf_1: cut_1) {
      for (auto leaf_2: cut_2) {
        if (leaf_1 == leaf_2) return_vector.push_back(leaf_1);
      }
    }
  }

  bool check_cut_legality(uint32_t i1_c0, uint32_t i1_c1, uint32_t i2_c0, uint32_t i2_c1,
          uint32_t cut_i1_1, uint32_t cut_i1_2, uint32_t cut_i2_1, uint32_t cut_i2_2) {

    auto cut_1_1 = cuts.cuts(i1_c0)[cut_i1_1];
    auto cut_1_2 = cuts.cuts(i1_c1)[cut_i1_2];
    auto cut_2_1 = cuts.cuts(i2_c0)[cut_i2_1];
    auto cut_2_2 = cuts.cuts(i2_c1)[cut_i2_2];
    uint32_t cut_size_1_1 = cut_1_1.size();
    uint32_t cut_size_1_2 = cut_1_2.size();
    uint32_t cut_size_2_1 = cut_2_1.size();
    uint32_t cut_size_2_2 = cut_2_2.size();

    // 4 pins to each LUT
    if ((cut_size_1_1 > 4)||(cut_size_1_2 > 4)||(cut_size_2_1 > 4)||(cut_size_2_2 > 4))
      return false;

    uint32_t nmatch_1_1_1_2 = number_of_match_in_cuts(cut_1_1, cut_1_2);
    uint32_t nmatch_1_1_2_1 = number_of_match_in_cuts(cut_1_1, cut_2_1);
    uint32_t nmatch_1_1_2_2 = number_of_match_in_cuts(cut_1_1, cut_2_2);
    uint32_t nmatch_1_2_2_1 = number_of_match_in_cuts(cut_1_2, cut_2_1);
    uint32_t nmatch_1_2_2_2 = number_of_match_in_cuts(cut_1_2, cut_2_2);
    uint32_t nmatch_2_1_2_2 = number_of_match_in_cuts(cut_2_1, cut_2_2);

    // 5 unique pins to top/bottom half of LUT
    if ((cut_size_1_1 + cut_size_1_2 - nmatch_1_1_1_2) > 5){
      return false;
    }
    if ((cut_size_2_1 + cut_size_2_2 - nmatch_2_1_2_2) > 5) {
      return false;
    }

    // how many a/b connections does top/bottom half of LUT need?
    uint32_t abc_1 = std::max(0, int(3-int(4-cut_size_1_1)));
    uint32_t abc_2 = std::max(0, int(3-int(4-cut_size_1_2)));
    uint32_t abd_1 = std::max(0, int(3-int(4-cut_size_2_1)));
    uint32_t abd_2 = std::max(0, int(3-int(4-cut_size_2_2)));

    if (int(abc_1 + abc_2 - nmatch_1_1_1_2) > 3) {
      return false;
    }
    if (int(abd_1 + abd_2 - nmatch_2_1_2_2) > 3) {
      return false;
    }
    
    // TODO: missing a/b shared pins condition
    return true;
  }

  ///////////////////////////////////////////////////////////////
  // LUT Mapping 
  ///////////////////////////////////////////////////////////////
  template<bool ELA>
  void compute_mapping( bool delay )
  {
    uint32_t time = 0;
    node<Ntk> slowest_leaf;
    for ( auto const& n : top_order )
    {
      if ( ntk.is_constant( n ) || ntk.is_pi( n ))
        continue;

      if (is_a_carry_node(n)) {
        time = 0;
        slowest_leaf = n;
        if (ps.verbose && ps.verbosity > 3) std::cout << "Node* " << ntk.node_to_index(n) << ":";
        for ( auto leaf : carry_cut_list[ntk.node_to_index(n)] )
        {
          //time = std::max( time, delays[leaf] );
          if (delays[leaf] > time) {
            time = delays[leaf];
            slowest_leaf = leaf;
          }
          if (ps.verbose && ps.verbosity > 3)std::cout << leaf << "(" << delays[leaf] << ") ";
        }
        if (carry_driver_nodes[n] == slowest_leaf)
          time += CARRY_DELAY;
        else //if (is_a_carry_node (slowest_leaf))
          time += LUT_ADDER_DELAY;

        map_refs[ntk.node_to_index(n)] = 0;
        delays[ntk.node_to_index(n)] = time;
        if (ps.verbose && ps.verbosity > 3)std::cout << "-> " << slowest_leaf << "(" << delays[slowest_leaf] << ")\n";
        continue;
      }  
      compute_best_cut<ELA>( ntk.node_to_index( n ), delay );
    }
    set_mapping_refs<ELA>();
  }

  template<bool ELA>
  void set_carry_mapping_refs()
  {

    /* compute current delay and update mapping refs */
    delay = 0;
    delay_lut = 0;
    ntk.foreach_po( [this]( auto s ) {
      const auto index = ntk.node_to_index( ntk.get_node( s ) );
      delay = std::max( delay, delays[index] );
      if (is_a_carry_node(ntk.get_node(s))) {
        map_refs[index] = 0;
        //delay_lut = std::max( delay_lut, delays[index]);
      }
    } );

    // Nodes driving the carry need to be mapped
    // Unless they are a carry itself or an input
    for (auto const& n : top_order) {
      if (!is_a_carry_node(n)) {
        delay_lut = std::max( delay_lut, delays[ntk.node_to_index(n)]);
      } 
    }
  }

  template<bool ELA>
  void set_mapping_refs()
  {

    const auto coef = 1.0f / ( 1.0f + ( iteration + 1 ) * ( iteration + 1 ) );

    /* compute current delay and update mapping refs */
    delay = 0;
    delay_lut = 0;
    ntk.foreach_po( [this]( auto s ) {
      const auto index = ntk.node_to_index( ntk.get_node( s ) );
      delay = std::max( delay, delays[index] );
      // only increase the map_ref if the PO is not fed by carry chain
      if (!is_a_carry_node(ntk.get_node(s))) {
        delay_lut = std::max( delay_lut, delays[index]);
        if constexpr ( !ELA ) {
          map_refs[index]++;
        }
      }
    } );

    // Nodes driving the carry need to be mapped
    // Unless they are a carry itself or an input
    for (auto const& n : top_order) {
      if (is_a_carry_node(n)) {
        for (auto const cut_index: carry_cut_list[ntk.node_to_index(n)]) {
          if (!is_a_carry_node(ntk.index_to_node(cut_index)) && \
             !ntk.is_pi(ntk.index_to_node(cut_index))) {
            map_refs[cut_index]++;
          }
        }
      } else {
        delay_lut = std::max( delay_lut, delays[ntk.node_to_index(n)]);
      }
    }
  
    /* compute current area and update mapping refs */
    area = 0;
    for ( auto it = top_order.rbegin(); it != top_order.rend(); ++it )
    {
      /* skip constants and PIs (TODO: stop earlier) */
      if ( ntk.is_constant( *it ) || ntk.is_pi( *it ) )
        continue;

      const auto index = ntk.node_to_index( *it );
      if ( map_refs[index] == 0 )
        continue;

      if constexpr ( !ELA )
      {
        for ( auto leaf : cuts.cuts( index )[0] )
        {
          if (!is_a_carry_node(leaf))
            map_refs[leaf]++;
        }
      }
      area++;
    }

    /* blend flow referenes */
    for ( auto i = 0u; i < ntk.size(); ++i )
    {
      flow_refs[i] = coef * flow_refs[i] + ( 1.0f - coef ) * std::max<float>( 1.0f, map_refs[i] );
    }
    
    ++iteration;
  }

  std::pair<float, uint32_t> cut_flow( cut_t const& cut)
  {
    uint32_t time{0u};
    float flow{0.0f};

    for ( auto leaf : cut )
    {
      time = std::max( time, delays[leaf] );
      // Penalty for choosing a cut that is driven by carry node (incurs additional cost)
      //if (is_a_carry_node(leaf)) {
      //  future_time = std::max( future_time, delays[leaf]+LUT_ADDER_DELAY );
      //  std::cout << "USING THIS\n";
      //  flow += (flows[leaf]*4);
      //} else {
        flow += flows[leaf];
      //}
    }

    return {flow + cut_area( cut ), time+LUT_DELAY};
  }

  uint32_t cut_area( cut_t const& cut ) const
  {
    return static_cast<uint32_t>( cut->data.cost );
  }


  /* reference cut:
   *   adds cut to current mapping and recursively adds best cuts of leaf
   *   nodes, if they are not part of the current mapping.
   */
  uint32_t cut_ref( cut_t const& cut )
  {
    uint32_t count = cut_area( cut );
    for ( auto leaf : cut )
    {
      if ( ntk.is_constant( ntk.index_to_node( leaf ) ) || ntk.is_pi( ntk.index_to_node( leaf ) ) )
        continue;

      if ( map_refs[leaf]++ == 0 )
      {
        count += cut_ref( cuts.cuts( leaf )[0] );
      }
    }
    return count;
  }

  /* dereference cut:
   *   removes cut from current mapping and recursively removes best cuts of
   *   leaf nodes, if they are part of the current mapping.
   *   (this is the inverse operation to cut_ref)
   */
  uint32_t cut_deref( cut_t const& cut )
  {
    uint32_t count = cut_area( cut );
    for ( auto leaf : cut )
    {
      if ( ntk.is_constant( ntk.index_to_node( leaf ) ) || ntk.is_pi( ntk.index_to_node( leaf ) ) )
        continue;

      if ( --map_refs[leaf] == 0 )
      {
        count += cut_deref( cuts.cuts( leaf ).best() );
      }
    }
    return count;
  }

  /* reference cut (special version):
   *   this special version of cut_ref does two additional things:
   *   1. it stops recursing if it has found `limit` cuts
   *   2. it remembers all cuts for which the reference count increases in the
   *      vector `tmp_area`.
   */
  uint32_t cut_ref_limit_save( cut_t const& cut, uint32_t limit )
  {
    uint32_t count = cut_area( cut );
    if ( limit == 0 )
      return count;

    for ( auto leaf : cut )
    {
      if ( ntk.is_constant( ntk.index_to_node( leaf ) ) || ntk.is_pi( ntk.index_to_node( leaf ) ) )
        continue;

      tmp_area.push_back( leaf );
      if ( map_refs[leaf]++ == 0 )
      {
        count += cut_ref_limit_save( cuts.cuts( leaf ).best(), limit - 1 );
      }
    }
    return count;
  }

  /* estimates the cost of adding this cut to the mapping:
   *   This algorithm references cuts recursively to estimate how many cuts
   *   would be needed to add to the mapping if `cut` were to be added.  It
   *   temporarily modifies the reference counters but reverts them eventually.
   */
  uint32_t cut_area_estimation( cut_t const& cut )
  {
    tmp_area.clear();
    const auto count = cut_ref_limit_save( cut, 8 );
    for ( auto const& n : tmp_area )
    {
      map_refs[n]--;
    }
    return count;
  }
 
  template<bool ELA>
  void compute_best_cut( uint32_t index, bool delay )
  {
    constexpr auto mf_eps{0.005f};

    float flow;
    uint32_t time{0};
    int32_t best_cut{-1};
    float best_flow{std::numeric_limits<float>::max()};
    uint32_t best_time{std::numeric_limits<uint32_t>::max()};
    int32_t cut_index{-1};

    if constexpr ( ELA )
    {
      if ( map_refs[index] > 0 )
      {
        cut_deref( cuts.cuts( index )[0] );
      }
    }

    for ( auto* cut : cuts.cuts( index ) )
    {
      ++cut_index;
      if ( cut->size() == 1 )
        continue;

      if constexpr ( ELA )
      {
        flow = static_cast<float>( cut_area_estimation( *cut ) );
      }
      else
      {
        std::tie( flow, time ) = cut_flow( *cut );
      }
      if ( delay ) {
        if ( best_cut == -1 ||  best_time > time ||  ( best_time == time && best_flow > flow - mf_eps ) )
        {
          best_cut = cut_index;
          best_flow = flow;
          best_time = time;
        }
      } else {
        if ( best_cut == -1 || best_flow > flow + mf_eps || ( best_flow > flow - mf_eps && best_time > time ) )
        {
          best_cut = cut_index;
          best_flow = flow;
          best_time = time;
        }
      }
    }

    if constexpr ( ELA )
    {
      if ( map_refs[index] > 0 )
      {
        cut_ref( cuts.cuts( index )[best_cut] );
      }
    }
    else
    {
      map_refs[index] = 0;
    }
    if constexpr ( ELA )
    {
      best_time = cut_flow( cuts.cuts( index )[best_cut] ).second;
    }
    delays[index] = best_time;
    flows[index] = best_flow / flow_refs[index];
    
    if ( best_cut != 0 )
    {
      cuts.cuts( index ).update_best( best_cut );
    }
    if (ps.verbose && ps.verbosity > 3)std::cout << "Node " << index << ": ";

    auto leaf_of_best_cut = 0;
    for (auto leaf: cuts.cuts(index)[0]) {
      if (ps.verbose && ps.verbosity > 3)std::cout << leaf << "(" << delays[leaf] << ") ";
      if (delays[leaf_of_best_cut] < delays[leaf]) {
        leaf_of_best_cut = leaf;
      }
    }
    if (ps.verbose && ps.verbosity > 3) {
      std::cout << "-> " << leaf_of_best_cut << "(" << delays[leaf_of_best_cut] << ")" ;
      std::cout << " : " << delays[index] <<"\n";
    }
  }

  ///////////////////////////////////////////////////////////////
  // Derive Mapping 
  ///////////////////////////////////////////////////////////////

  void derive_carry_mapping() {

    // For each path, add nodes to KLUT mapping
    // Set mapping index to the carry it's driven by
    for ( auto const path : carry_paths) {

      for (uint32_t i = 1; i < path.size(); i++) {

        auto const n = path[i];
        auto const ncarry = path[i-1];

        if ( ntk.is_constant( n ) || ntk.is_pi( n ) || ntk.is_cell_root( n ) )  
          continue;

        const auto index = ntk.node_to_index( n );

        std::vector<node<Ntk>> nodes;

        ntk.add_to_carry_mapping (n, ntk.node_to_index(ncarry)); 
        for ( auto c: carry_cut_list[ntk.node_to_index(n)]) {
          nodes.push_back(c);
        }
        ntk.add_to_mapping( n, nodes.begin(), nodes.end() );
      }
    }
  }

  void derive_mapping()
  {
    ntk.clear_mapping();

    derive_carry_mapping();

    for ( auto const& n : top_order )
    {
      if ( ntk.is_constant( n ) || ntk.is_pi( n ) || ntk.is_cell_root( n ) ) 
        continue;

      const auto index = ntk.node_to_index( n );

      std::vector<node<Ntk>> nodes;

      if ( map_refs[index] == 0 ) 
        continue;
      for ( auto const& l : cuts.cuts( index ).best() ) {
        nodes.push_back( ntk.index_to_node( l ) );
      }
      ntk.add_to_mapping( n, nodes.begin(), nodes.end() );

      // This is translated into LUT functionality
      if constexpr ( StoreFunction )
      {
        if ( !is_a_carry_node(n) ) {
          ntk.set_cell_function( n, cuts.truth_table( cuts.cuts( index ).best() ) );
        }
      }
    }
  }

  ///////////////////////////////////////////////////////////////
  // Supporting Functions 
  ///////////////////////////////////////////////////////////////

  bool is_a_carry_node(node<Ntk> n) {
    if (carry_nodes[n] != 0) return true;
    return false;
  }

  bool is_in_carry_lut(node<Ntk> n) {
    if (!carry_cut_list[ntk.node_to_index(n)].empty()) return true;
    return false;
  }

  void find_num_nodes(node<Ntk> n, cut_t const& bestcut, uint32_t* num_nodes,
      uint32_t depth, uint32_t* max_depth, bool printtt) {

    ntk.set_visited(n,1);

    for (uint32_t i = 0; i < ntk.fanin_size(n); i++) {
      bool found = false;
      auto c = ntk.get_children(n,i);
      for ( auto const& l : bestcut ) {
        if (l == ntk.node_to_index(c)) {
          if (*max_depth < depth)
            *max_depth = depth;
          found = true;
        } 
      } 
      if (!found && !ntk.is_pi(c) && !ntk.is_constant(c) && !(ntk.visited(c)>0)) {
        (*num_nodes)++;
        find_num_nodes(c, bestcut, num_nodes, depth+1, max_depth, printtt);
      }
    }
  }


  ///////////////////////////////////////////////////////////////
  // Print Functions
  ///////////////////////////////////////////////////////////////
  
  void print_cut_enumerations() {
    ntk.foreach_node( [this]( auto n, auto ) {
      std::cout << n << ":\n"; 
      for ( auto* cut : cuts.cuts( ntk.node_to_index(n)) ) {
        std::cout << "\t{ ";
        uint32_t flow_total = 0;
        for ( auto leaf: *cut) {
          std::cout << leaf << "(" << delays[leaf] << ") ";
          flow_total += flow_refs[leaf];
        }
        std::cout << "} (" << flow_total << ")\n";
      }
      std::cout << "\n";
    });
  }

  void print_inverter_from_path (void) {
    ntk.foreach_node( [this]( auto n, auto ) {
      std::cout << n << ":"; 
      for (uint32_t i = 0; i < ntk.fanin_size(n); i++) {
        node<Ntk> child_node = ntk.get_children(n,i);  
        std::cout << "  " << child_node << "(" << ntk.is_complemented_children(n,i) << ") ";
      }
      std::cout << "\n";
    });
  }

  void print_carry_nodes() {
    std::cout << "Carry Nodes: ";
    for (auto n_carry: carry_nodes) {
      std::cout << n_carry << ",";
    }
    std::cout << "\n";
  }

  void print_carry_paths() {
    std::cout << "Carry Paths\n";
    uint32_t count = 0;

    for (auto carry_path: carry_paths) {
      std::cout << count << ": ";
      for (auto const node: carry_path) {
        std::cout << node << ",";
      }
      count++;
      std::cout << "\n";
    }
  }

  void print_path (std::vector<node<Ntk>> path) {
    std::cout << "Path: ";
    for (auto const node: path) {
      if (is_in_carry_lut(node)) std::cout << "*";
      std::cout << node << ",";
    }
    std::cout << "\n";
  }

  void print_carry_lut_nodes() {
    std::cout << "Carry LUT Nodes\n";
    for (int i = 0; i < carry_cut_list.size(); i++) {
      std::cout << i << ": ";
      for (auto const& cut: carry_cut_list[i]) {
        std::cout << cut << ",";
      }
      std::cout << "\n";
    }
  }

  void print_mig_path_and_depth_helper (uint32_t n, uint32_t& path_count) {
    if (path_count > 1000000) {
      return;
    } else if (ntk.is_pi(n) || ntk.is_constant(n)) {
      path_count++; 
      return;
    } else {
      ntk.foreach_fanin( n, [&]( auto const& f ) {
        auto n_child = ntk.get_node(f);
        if (path_count < 1000000) {
          print_mig_path_and_depth_helper (n_child, path_count);
        } else return;
      });
    }
  }

  void print_LUT_path_and_depth_helper (uint32_t n, uint32_t& path_count) {
    uint32_t worst_delay = 0;
    if (path_count > 1000000) {
      return;
    } else if (ntk.is_pi(n) || ntk.is_constant(n)) {
      path_count++; 
      return;
    } else {
      auto index = ntk.node_to_index(n); 
      for (auto leaf: cuts.cuts(index)[0]) {
        if (delays[leaf] > worst_delay) {
          worst_delay = delays[leaf];
        }
      }
      for (auto leaf: cuts.cuts(index)[0]) {
        if (delays[leaf] == worst_delay) {
          print_LUT_path_and_depth_helper (leaf, path_count);
        }
      }
    }
  }

  uint32_t num_worst_paths() {
    uint32_t worst_critical_delay = 0;
    ntk.foreach_po( [&]( auto s ) {
      const auto n = ntk.get_node( s );
      if (worst_critical_delay < delays[n]) {
        worst_critical_delay = delays[n];
      }
    });

    uint32_t critical_path_count = 0;
    uint32_t path_count = 0;
    uint32_t total_path_count = 0;
    uint32_t total_critical_path_count = 0;

    ntk.foreach_po( [&]( auto s ) {
      critical_path_count = 0;
      path_count = 0;
      const auto n = ntk.get_node( s );
      print_mig_path_and_depth_helper (n, path_count);
      total_path_count += path_count;
      if (delays[n] == worst_critical_delay) {
        print_LUT_path_and_depth_helper (n, critical_path_count);
        std::cout << "PO " << n << "(" << delays[n] << ") has " << critical_path_count << " critical paths.\n"; 
        total_critical_path_count += critical_path_count;
      }
    });
    std::cout << "Total critical path is " << total_critical_path_count << " out of " <<  total_path_count << "\n";
    return total_critical_path_count;
  }

  void print_critical_path_helper (uint32_t index, uint32_t* lut_count, uint32_t* carry_count) {

    if (is_a_carry_node(ntk.index_to_node(index))) {
      std::cout << " -> " << index << "*(" << delays[index] << "): ";
      (*carry_count)++;
    } else {
      std::cout << " -> " << index << "(" << delays[index] << "): ";
      (*lut_count)++;
    }

    if (ntk.is_pi(ntk.index_to_node(index)) || ntk.is_constant(ntk.index_to_node(index)))
      return;
   
    uint32_t worst_delay = 0;
    uint32_t critical_index = 0;
    if (is_a_carry_node(ntk.index_to_node(index))) {
      for ( auto leaf : carry_cut_list[index] ) {
        std::cout << leaf << " ";
        if (delays[leaf] > worst_delay) {
          critical_index = leaf;
          worst_delay = delays[leaf];
        }
      }
    } else {
      for (auto leaf: cuts.cuts(index)[0]) {
        std::cout << leaf << " ";
        if (delays[leaf] > worst_delay) {
          critical_index = leaf;
          worst_delay = delays[leaf];
        }
      }
    }
    std::cout << "\n";
    print_critical_path_helper (critical_index, lut_count, carry_count);

  }

  // Print the critical MIG path
  void print_critical_path () {

    std::cout << "Critical path:\n";
  
    uint32_t worst_delay = 0;
    uint32_t critical_index = 0;
    uint32_t lut_count = 0;
    uint32_t carry_count = 0;
    ntk.foreach_po( [&]( auto s ) {
      const auto index = ntk.node_to_index( ntk.get_node( s ) );
      if (delays[index] > worst_delay) {
        critical_index = index;
        worst_delay = delays[index];
      }
    });
    print_critical_path_helper (critical_index, &lut_count, &carry_count);
    std::cout << "total: " << lut_count << " " << carry_count << "\n";
  }


  void print_state()
  {
    if (ps.verbosity > 1) {
      for ( auto i = 0u; i < ntk.size(); ++i )
      {
        std::cout << fmt::format( "*** Obj = {:>3} (node = {:>3})  FlowRefs = {:5.2f}  MapRefs = {:>2}  Flow = {:5.2f}  Delay = {:>3} Carry = {}\n", \
          i, ntk.index_to_node( i ), flow_refs[i], map_refs[i], flows[i], delays[i], carry_nodes[i]);
      }
    }
    std::cout << fmt::format( "Level = {}  Area = {}\n", delay/LUT_DELAY, area );
  }

private:
  Ntk& ntk;
  carry_lut_mapping_params const& ps;
  carry_lut_mapping_stats& st;

  uint32_t iteration{0}; /* current mapping iteration */
  uint32_t delay{0};     /* current delay of the mapping */
  uint32_t area{0};      /* current area of the mapping */
  uint32_t delay_lut{0};
  //bool ela{false};       /* compute exact area */

  std::vector<node<Ntk>> top_order;
  std::vector<float> flow_refs;
  std::vector<uint32_t> map_refs;
  std::vector<float> flows;
  std::vector<uint32_t> delays;
  std::vector<uint32_t> paths;
  std::vector<std::vector<uint32_t>> carry_cut_list;
  std::vector<node<Ntk>> carry_nodes;
  std::vector<node<Ntk>> carry_driver_nodes;
  std::vector<std::vector<node<Ntk>>> carry_paths;
  network_cuts_t cuts;

  std::vector<uint32_t> tmp_area; /* temporary vector to compute exact area */

  // Contains the list of nodes to be mapped to carry 
};

}; /* namespace detail */


/*! \brief LUT mapping.
 *
 * This function implements a LUT mapping algorithm.  It is controlled by two
 * template arguments `StoreFunction` (defaulted to `true`) and `CutData`
 * (defaulted to `cut_enumeration_mf_cut`).  The first argument `StoreFunction`
 * controls whether the LUT function is stored in the mapping.  In that case
 * truth tables are computed during cut enumeration, which requires more
 * runtime.  The second argument is simuilar to the `CutData` argument in
 * `cut_enumeration`, which can specialize the cost function to select priority
 * cuts and store additional data.  For LUT mapping using this function the
 * type passed as `CutData` must implement the following three fields:
 *
 * - `uint32_t delay`
 * - `float flow`
 * - `float costs`
 * 
 * See `include/mockturtle/algorithms/cut_enumeration/mf_cut.hpp` for one
 * example of a CutData type that implements the cost function that is used in
 * the LUT mapper `&mf` in ABC.
 *
 * **Required network functions:**
 * - `size`
 * - `is_pi`
 * - `is_constant`
 * - `node_to_index`
 * - `index_to_node`
 * - `get_node`
 * - `foreach_po`
 * - `foreach_node`
 * - `fanout_size`
 * - `clear_mapping`
 * - `add_to_mapping`
 * - `set_lut_funtion` (if `StoreFunction` is true)
 *
   \verbatim embed:rst

   .. note::

      The implementation of this algorithm was heavily inspired but the LUT
      mapping command ``&mf`` in ABC.
   \endverbatim
 */
template<class Ntk, bool StoreFunction = false, typename CutData = cut_enumeration_mf_cut>
void carry_lut_mapping( Ntk& ntk, carry_lut_mapping_params const& ps = {}, carry_lut_mapping_stats* pst = nullptr )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_size_v<Ntk>, "Ntk does not implement the size method" );
  static_assert( has_is_pi_v<Ntk>, "Ntk does not implement the is_pi method" );
  static_assert( has_is_constant_v<Ntk>, "Ntk does not implement the is_constant method" );
  static_assert( has_node_to_index_v<Ntk>, "Ntk does not implement the node_to_index method" );
  static_assert( has_index_to_node_v<Ntk>, "Ntk does not implement the index_to_node method" );
  static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
  static_assert( has_foreach_po_v<Ntk>, "Ntk does not implement the foreach_po method" );
  static_assert( has_foreach_node_v<Ntk>, "Ntk does not implement the foreach_node method" );
  static_assert( has_fanout_size_v<Ntk>, "Ntk does not implement the fanout_size method" );
  static_assert( has_clear_mapping_v<Ntk>, "Ntk does not implement the clear_mapping method" );
  static_assert( has_add_to_mapping_v<Ntk>, "Ntk does not implement the add_to_mapping method" );
  static_assert( !StoreFunction || has_set_cell_function_v<Ntk>, "Ntk does not implement the set_cell_function method" );

  carry_lut_mapping_stats st;
  detail::carry_lut_mapping_impl<Ntk, StoreFunction, CutData> p( ntk, ps, st );
  p.run();
  if ( ps.verbose )
  {
    st.report();
  }

  if ( pst )
  {
    *pst = st;
  }
}

}
