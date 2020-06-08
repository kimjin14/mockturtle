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
  uint32_t rounds{2u};

  /*! \brief Number of rounds for exact area optimization. */
  uint32_t rounds_ela{1u};

  /*! \brief Be verbose. */
  bool verbose{true};

  /*! \brief Verbosity Level. >3 means print everything*/
  uint32_t verbosity = 4;

  /* Determines whether to print carry node combined LUT fcn. SHOULD NOT BE USED ANYMORE! */
  bool carry_lut_combined{false};  

  /* Map to carry. */
  bool carry_mapping{true};

  /* Selects to use Xilinx architecture instead of Intel. */
  bool xilinx_arch{true};

  /* Number of rounds for carry chain synthesis. */
  uint32_t max_rounds_carry{300};  
 
  /* Cost function to be used. */ 
  int cost{1};

  /* 5-LUT sharing. */
  bool lut_sharing{true};

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
        carry_cut_index_list( ntk.size() ),
        carry_nodes( ntk.size(), 0 ),
        carry_driver_nodes( ntk.size(), 0 ),
        mapped_to_5LUTa( ntk.size(), 0 ),
        mapped_to_5LUTb( ntk.size(), 0 ),
        mapped_to_5LUT( ntk.size(), false),
        mapped_to_5LUT_complemented( ntk.size(), false ),
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
      print_state();
    }

    /* map to carry chain if true */
    if (ps.carry_mapping)
    {
      carry_chain_mapping();
      print_state();
    }

    while ( iteration < ps.rounds + ps.rounds_ela)
    {
      compute_mapping<true>(false);
      print_state();
    }

    //carry_chain_mapping_test();
    if (ps.carry_mapping) update_delay();
    print_critical_path();
    print_state();

    derive_mapping();

    std::cout << "PRINTING INVERTER\n";
    uint32_t counting = 0;
    for (auto eachnode: mapped_to_5LUT_complemented) {
      if (eachnode != 0)
        std::cout << counting << ": " << eachnode << "\n";
      counting++;
    }
  }

private:

  void carry_chain_mapping_test () {

      // Determine how many rounds of carry chain mapping
      uint32_t rounds_carry_mapping = 10;
      rounds_carry_mapping = num_worst_paths();
      if (rounds_carry_mapping > ps.max_rounds_carry)
        rounds_carry_mapping = ps.max_rounds_carry;
  
      // Iteratively map paths to carry chain
      for (uint32_t i = 0; i < rounds_carry_mapping; i++) { 

        if (ps.verbose) std::cout << "Carry Mapping Iteration " \
          << i << " targetting delay of " << delay << ".\n";

        select_path_and_map();

        update_delay();   
        print_state();
      }
      print_state();
      std::cout << "There were " << carry_paths.size() << " paths placed on carry chain\n"; 
  }


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

        uint32_t delay_offset = 2 * LUT_DELAY;
        if (ps.verbose) std::cout << "Carry Mapping Iteration " \
          << i << " targetting delay of " << delay << " with offset " << delay_offset << ".\n";

        // Select nodes to be placed on carry 
        std::vector<node<Ntk>> path_for_carry_chain;
        if (!path_selection(path_for_carry_chain, delay_offset))
          break;

        // Compute carry LUT mapping 
        if (ps.xilinx_arch)
          xilinx_compute_carry_mapping(path_for_carry_chain);
        else
          compute_carry_mapping(path_for_carry_chain);

        update_delay();   
        path_for_carry_chain.clear();
        print_state();
      }

      set_carry_mapping_refs();
      remove_inverter_for_carry_mapping();
      print_state();
      std::cout << "There were " << carry_paths.size() << " paths placed on carry chain\n"; 
  }

  ///////////////////////////////////////////////////////////////
  // Selecting Path to be Placed on Carry Chain and Mapping
  ///////////////////////////////////////////////////////////////
  uint32_t count_path_to_node (std::vector<node<Ntk>>& path_for_carry_chain, 
      uint32_t index, uint32_t source_index, uint32_t dest_index, uint32_t& node_length) {

    auto const& n = ntk.index_to_node(source_index);

    if (source_index == dest_index) {
      return 1;
    }

    if (is_a_carry_node(ntk.index_to_node(index))) {
      for ( auto leaf : carry_cut_list[index] )
        if (source_index == leaf) return 0;
    } else {
      for ( auto leaf : cuts.cuts(index)[0] )
        if (source_index == leaf) return 0;
    }
 
    uint32_t total = 0;

    for (uint32_t i = 0; i < ntk.fanin_size(n); i++) {
      auto const leaf = ntk.get_children(n,i);
      auto const leaf_index = ntk.node_to_index(leaf);
      if (!is_a_carry_node(leaf) && !ntk.is_constant(leaf))
        total += count_path_to_node(path_for_carry_chain, index, leaf_index, dest_index, node_length );
    }
    return total;
  }

  bool add_node_and_cut_to_LUT  (std::vector<node<Ntk>>& path_for_carry_chain, 
      uint32_t index, uint32_t source_index, uint32_t dest_index, uint32_t& node_length) {

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
        found = add_node_and_cut_to_LUT (path_for_carry_chain, index, leaf_index, dest_index, node_length );
      if (found) {
        node_length++;
        if (ps.verbose && ps.verbosity > 2) std::cout << "adding " << leaf << "\n";

        // Add node to path
        path_for_carry_chain.push_back(leaf);

        // Update cut list for the carry node
        std::cout << "\tCut list: ";
        uint32_t max_delay = 0;
        for (auto cut_leaf: cuts.cuts(index)[0]) {
          if (cut_leaf != dest_index) {
            std::cout << cut_leaf << " ";
            carry_cut_list[source_index].push_back(cut_leaf);
            if (delays[cut_leaf] > max_delay) max_delay = delays[cut_leaf];
          }
        } 
        carry_cut_list[source_index].push_back(leaf_index);
        std::cout << leaf_index << "\n";

        // Update delay
        std::cout << "\tUpdate delay for node " << source_index;
        if (map_refs[source_index] <= 0) std::cout << "*";  
        std::cout << ": " << delays[source_index];
        delays[source_index] = std::max(max_delay + LUT_ADDER_DELAY, delays[leaf_index] + CARRY_DELAY);
        std::cout << " -> " << delays[source_index] << "\n";

        // Return that you can map this node to carry
        return true;
      }
    }
    return found;
  }

  //cut_t find_children_cut (uint32_t index, uint32_t c_index1) {
  //  
  //
  //}

  // Determine which leaf should be considered for carry chain
  uint32_t select_leaf( auto cut_list ) {

    uint32_t max_leaf = 0;
    uint32_t max_delay = 0;
    float min_fanout{std::numeric_limits<float>::max()};

    // Find worst delay
    for ( auto leaf : cut_list ) {
      if (ps.verbose && ps.verbosity > 2) std::cout << "\t" << leaf << ":" << delays[leaf] << "," << flow_refs[leaf] << "," << ntk.fanout_size(leaf) <<"\n"; 
      if (max_delay < delays[leaf] || (max_delay == delays[leaf] && (ntk.fanout_size(leaf) < min_fanout))) {
        max_leaf = leaf;
        max_delay = delays[leaf];
        min_fanout = flows[leaf];
      }
    }

    return max_leaf;
  }

  // Keep going deeper to its fanin until depth is found
  bool find_LUT_for_reduction (std::vector<node<Ntk>>& path_for_carry_chain, uint32_t index, uint32_t length) { 
    
    if (ps.verbose && ps.verbosity > 2) std::cout << "Index " << index << "(" << delays[index] << "):";
    //if (delays[index] == 0 || ntk.is_pi(index)) {
    if (delays[index] == 0 && length > 4) {
      if (ps.verbose && ps.verbosity > 2) std::cout << "\n";
      return true;
    }

    bool deepest = false;

    uint32_t max_leaf = index;

    if (is_a_carry_node(ntk.index_to_node(index))) {
      if (ps.verbose && ps.verbosity > 2) std::cout << " carry\n"; 
      max_leaf = select_leaf(carry_cut_list[index]);
    } else {
      if (ps.verbose && ps.verbosity > 2) std::cout << " lut\n"; 
      max_leaf = select_leaf(cuts.cuts(index)[0]);
    }

    if (max_leaf != index) {
      deepest = find_LUT_for_reduction (path_for_carry_chain, max_leaf, length+1);
      if (deepest)  {
        uint32_t node_length = 0;
        std::cout << "count is " << count_path_to_node (path_for_carry_chain, index, index, max_leaf, node_length) << "\n";
        if(add_node_and_cut_to_LUT(path_for_carry_chain, index, index, max_leaf, node_length)) {
          //if (node_length == 1) {
          //  std::cout << "\t\tShould add cut here\n";
          //  path_for_carry_chain.push_back(index);
          //  for (auto cut_leaf: cuts.cuts(index)[0]) {
          //    if (cut_leaf != max_leaf) carry_cut_list[index].push_back(cut_leaf);
          //  } carry_cut_list[index].push_back(max_leaf);
          //}
          //if (ps.verbose && ps.verbosity > 2) std::cout << "\t\tAdded path of length " << node_length << "\n";
        } else { //cannot put this in carry
          if (ps.verbose && ps.verbosity > 2) std::cout << "\t\tCannot add to path\n";
          return false;
        } 
      }
    } 

    return deepest;
  }

  void select_path_and_map () {

    bool reduce_possible = true;
    std::vector<node<Ntk>> path_for_carry_chain;
    uint32_t delay_offset = 2*LUT_DELAY;

    while (reduce_possible) {
  
      for (uint32_t curr_offset = 0; curr_offset <= delay_offset; curr_offset+=LUT_DELAY) {
        for ( auto it = top_order.rbegin(); it != top_order.rend(); ++it ) {
          const auto node = *it;
          const auto index = ntk.node_to_index(node);

          // Node with worst delay
          if (delays[index] >= delay-curr_offset&& delay != 0 && !ntk.is_pi(node)) {
            std::cout << "Found target index " << index << "(" << delays[index] << ")\n";

            // Try to place a path starting from this node
            if (find_LUT_for_reduction(path_for_carry_chain,index, 0)) {
              if (!is_a_carry_node(index)) {
                //std::cout << "adding " << index << "\n";
                path_for_carry_chain.push_back(index);
                //carry_nodes[index] += 1;
              }
              // Only place one path at a time
              break;
            } else {
              path_for_carry_chain.clear();
            }
          }
        }
        if (!path_for_carry_chain.empty()) break;
      }
      if (path_for_carry_chain.empty()) break;

      print_path(path_for_carry_chain);

      // Keep the path information in carry_paths
      carry_paths.push_back(path_for_carry_chain);

      // Update driver info for figuring out which node drives which carry
      // in case there are 2 nodes being driving a carry node
      for (uint32_t j = 1; j < path_for_carry_chain.size(); j++) {
        carry_nodes[path_for_carry_chain[j]] += 1;
        carry_driver_nodes[path_for_carry_chain[j]] = path_for_carry_chain[j-1];
      }

      // Clear path and decide whether to run again
      path_for_carry_chain.clear();
      reduce_possible = false;
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
        if (ps.verbose && ps.verbosity > 2) std::cout << "adding " << leaf << "\n";
        return true;
      }
    }
    return found;
  }

  // Keep going deeper to its fanin until depth is found
  bool find_deepest_LUT (std::vector<node<Ntk>>& path_for_carry_chain, uint32_t index, uint32_t length) { 
    
    if (ps.verbose && ps.verbosity > 2) std::cout << "Index " << index << "(" << delays[index] << "):";
    //if (delays[index] == 0 || ntk.is_pi(index)) {
    if (delays[index] == 0 && length > 6) {
      if (ps.verbose && ps.verbosity > 2) std::cout << "\n";
      return true;
    }

    bool deepest = false;

    uint32_t max_leaf = index;
    uint32_t max_delay = 0;
    float min_fanout{std::numeric_limits<float>::max()};

    if (is_a_carry_node(ntk.index_to_node(index))) {
      if (ps.verbose && ps.verbosity > 2) std::cout << " carry\n"; 
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
      if (ps.verbose && ps.verbosity > 2) std::cout << " lut\n"; 
      for ( auto leaf : cuts.cuts( index )[0] ){
        if (!is_a_carry_node(leaf)) { 
          if (ps.verbose && ps.verbosity > 2) std::cout << "\t" << leaf << ":" << delays[leaf] << "," << flow_refs[leaf] << "," << ntk.fanout_size(leaf) <<"\n"; 
          if (max_delay < delays[leaf] || (max_delay == delays[leaf] && (ntk.fanout_size(leaf) < min_fanout))) {
            max_leaf = leaf;
            max_delay = delays[leaf];
            min_fanout = ntk.fanout_size(leaf);
          }
        }
      }
    }

    if (max_leaf != index) {
      deepest = find_deepest_LUT(path_for_carry_chain, max_leaf, length+1);

      if (deepest)  {
        uint32_t node_length = 0;
        std::cout << "count is " << count_path_to_node (path_for_carry_chain, index, index, max_leaf, node_length) << "\n";
        ntk.clear_visited();
        if(add_find_longest_path_in_mapping(path_for_carry_chain, index, index, max_leaf)) {
          if (ps.verbose && ps.verbosity > 2) std::cout << "\t\tAdded path\n";
        } else { //cannot put this in carry
          return false;
        } 
      }
    } //else assert(0);

    return deepest;
  }
  
  bool path_selection (std::vector<node<Ntk>>& path_for_carry_chain, uint32_t delay_offset) {

    for (uint32_t curr_offset = 0; curr_offset <= delay_offset; curr_offset+=LUT_DELAY) {
      for ( auto it = top_order.rbegin(); it != top_order.rend(); ++it ) {
        const auto node = *it;
        const auto index = ntk.node_to_index(node);

        // Node with worst delay
        if (delays[index] >= delay-curr_offset && delay != 0 /*&& map_refs[index] > 0 */ && !ntk.is_pi(node)) {
          std::cout << "Found target index " << index << "(" << delays[index] << ")\n";

          // Try to place a path starting from this node
          if (find_deepest_LUT(path_for_carry_chain,index, 0)) {
            if (!is_a_carry_node(index)) {
              //std::cout << "adding " << index << "\n";
              path_for_carry_chain.push_back(index);
              //carry_nodes[index] += 1;
            }
            // Only place one path at a time
            break;
          } else {
            path_for_carry_chain.clear();
          }
        }
      }
      if (!path_for_carry_chain.empty()) break;
    }
    if (path_for_carry_chain.empty()) return false;

    print_path(path_for_carry_chain);

    // Keep the path information in carry_paths
    carry_paths.push_back(path_for_carry_chain);

    // Update driver info for figuring out which node drives which carry
    // in case there are 2 nodes being driving a carry node
    for (uint32_t j = 1; j < path_for_carry_chain.size(); j++) {
      carry_nodes[path_for_carry_chain[j]] += 1;
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

      // Remove inverters in its path
      remove_inverter();

      // Must update cut truth tables since inverters were moved
      // This must be done after all inverters are removed
      for (auto n: top_order) {
        uint32_t index = ntk.node_to_index(n);
        if (is_a_carry_node(n)) { 
          std::cout << "CARRY NODE " << n << ":\n";
          carry_truth_table (n, carry_driver_nodes[index]);
        }
      }

      check_inverter();
  }


  void carry_truth_table (node<Ntk> n, node<Ntk> nc) {

    auto i = ntk.node_to_index(n);
    auto ic = ntk.node_to_index(nc);

    assert (i != 0);

    // Get relevant children node for mapping
    node<Ntk> i_child[2] = {0};
    get_children_node (i_child, i, ic); 

    cut_t const& cut_1 = cuts.cuts(i_child[0])[carry_cut_index_list[i][0]];
    cut_t const& cut_2 = cuts.cuts(i_child[1])[carry_cut_index_list[i][1]];

    // Set LUT function
    // cut_1 and cut_2 with MIG carry node should be
    // turned into truth table 
    uint32_t unique_cuts[5] = {0};
    uint32_t n_total = 0;   

    n_total = carry_cut_list[i].size() - 1;
    uint32_t curr_i = 0;
    for (auto leaf: carry_cut_list[i]) {
      if (curr_i == (n_total)) break;
      unique_cuts[curr_i] = leaf;
      curr_i++;
    }
  
    bool child_complement[3] = {0};
    determine_child_complement (child_complement, n, ic, i_child[0], i_child[1]);
    kitty::dynamic_truth_table function = compute_carry_function(carry_cut_index_list[i][0], carry_cut_index_list[i][1], i_child[0], i_child[1], child_complement, unique_cuts, n_total, cut_1, cut_2);
    ntk.set_cell_function(n, function);

  }

  // Remove inverters in the carry path
  // 32 -> 193 -> !193 -> 249
  // inputs of 193 will be inverted and 193 won't be anymore
  void remove_inverter() {
    for (auto carry_path: carry_paths) {
      for (uint32_t carry_i = 1; carry_i < carry_path.size(); carry_i++) {

        auto carry_node = carry_path[carry_i];
        auto carry_child_node = carry_path[carry_i-1];
        if (carry_child_node == carry_driver_nodes[carry_i]) assert(0);
        
        bool complemented_carry_edge = false;
        for (uint32_t i = 0; i < ntk.fanin_size(carry_node); i++) {
          auto child_node = ntk.get_children(carry_node,i);  
          if (child_node == carry_child_node && ntk.is_complemented_children(carry_node,i)) {
            complemented_carry_edge = true;
          }
        }

        // If the carry path is complemented, go through all the nodes
        // and flip its children that match this node and the node itself
        if (complemented_carry_edge) {

          if (ps.verbose && ps.verbosity > 3) std::cout << "\tFor node " << carry_node << "\n";

          // flip children
          for (uint32_t i = 0; i < ntk.fanin_size(carry_node); i++) {
            node<Ntk> child_node = ntk.get_children(carry_node,i);  
            if (ps.verbose && ps.verbosity > 3) std::cout << "\t\tflipping child " << child_node << "\n";
            if (ps.verbose && ps.verbosity > 3) std::cout << "\t\t\tfrom " << ntk.is_complemented_children(carry_node,i); 

            ntk.flip_children(carry_node, i);

            if (ps.verbose && ps.verbosity > 3) std::cout << " to " << ntk.is_complemented_children(carry_node,i) << "\n"; 
          }
          
          ntk.foreach_node( [&]( auto n, auto ) {
            for (uint32_t i = 0; i < ntk.fanin_size(n); i++) {
              if (ntk.get_children(n,i) == carry_node) {\
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
        } // completed_carry_edge

        for (uint32_t i = 0; i < ntk.fanin_size(carry_node); i++) {
          auto child_node = ntk.get_children(carry_node,i);  
          if (!is_a_carry_node(child_node) && mapped_to_5LUTa[carry_node] != 0 && mapped_to_5LUTa[carry_node] == child_node) {
            mapped_to_5LUT_complemented[child_node] = ntk.is_complemented_children(carry_node,i);
            std::cout << "LUT a child " << child_node << "'s complemented is " << ntk.is_complemented_children(carry_node,i) << "\n";
          }
          if (!is_a_carry_node(child_node) && mapped_to_5LUTb[carry_node] != 0 && mapped_to_5LUTb[carry_node] == child_node) {
            mapped_to_5LUT_complemented[child_node] = ntk.is_complemented_children(carry_node,i);
            std::cout << "LUT b child " << child_node << "'s complemented is " << ntk.is_complemented_children(carry_node,i) << "\n";
          }
        }
      }
    }
  }

  ///////////////////////////////////////////////////////////////
  // Initialization
  ///////////////////////////////////////////////////////////////

  void init_carry_chain_mapping() {
    //do nothing for now
    ntk.foreach_node( [this]( auto n, auto ) {
      const auto index = ntk.node_to_index( n );
      map_refs[index] = 0;
    });
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

  void xilinx_compute_carry_mapping (std::vector<node<Ntk>> path_for_carry_chain)
  {

    node<Ntk> n_carryin, n_first, n_second;

    // This function should not be called with an empty carry path
    assert(!path_for_carry_chain.empty());

    // Map nodes in pairs (i+=2) since each ALM has 2 adders
    for (uint32_t i = 1; i < path_for_carry_chain.size(); i++) {

      // Assign carryin, 1st/2nd node
      n_carryin = path_for_carry_chain[i-1];
      n_first = path_for_carry_chain[i];

      if (ps.verbose && ps.verbosity > 2) {
        std::cout << "Carry LUT mapping iteration " << i << ": " << n_carryin << " " << n_first << "\n";
      }

      xilinx_compute_carry_LUT_mapping(n_first, n_carryin);
    }
  }

  // This function takes 2 index with its carry-in index and figures out 
  // if it cuts can be mapped to LUTs before carry node and the cost
  int xilinx_compute_carry_LUT_mapping(uint32_t i, uint32_t ic) {

    double cost = 0;
    double min_cost = std::numeric_limits<double>::max();
    int32_t max_i_1 = -1;
    int32_t max_i_2 = -1;
    constexpr auto eps{0.00005f};

    // Get relevant children node for mapping
    node<Ntk> i_child[2] = {0};
    get_children_node (i_child, i, ic); 

    uint32_t total_num_cuts = 0;
    uint32_t total_num_legal_cuts = 0;

    for (uint32_t cut_i_1 = 0; cut_i_1 < cuts.cuts(i_child[0]).size(); cut_i_1++) {
      for (uint32_t cut_i_2 = 0; cut_i_2 < cuts.cuts(i_child[1]).size(); cut_i_2++) {

        total_num_cuts++;

        // Array holding the nodes to the halfs of ALM as separate 4-LUT
        // [0][0] and [1][0] contains # of input
        uint32_t i_child_cuts[2][5] = {0}; 
        uint32_t abce0f0[5] = {0};
  
        ////////////////////////////////
        // Checking input requirements
        //////////////////////////////// 
        bool check_legality = true;
        if (check_xilinx_5lut_legality (i_child[0], i_child[1], cut_i_1, cut_i_2)) {
        } else check_legality = false;

        if (!check_legality) continue;

        total_num_legal_cuts++;

        ////////////////////////////////
        // Calculating cost of cut 
        //////////////////////////////// 

        cost = cost_of_5LUT_cuts(i, i_child[0], i_child[1], cut_i_1, cut_i_2);

        // The new cost should be greater than the max cost + floating-point noise 
        if (cost < min_cost - eps) {
          min_cost = cost;
          max_i_1 = cut_i_1;
          max_i_2 = cut_i_2;
        }
        if (ps.verbose && ps.verbosity > 2) {
          print_cut_cost (i, i_child[0], i_child[1], cut_i_1, cut_i_2, cost); 
        }
      }
    }
    if (ps.verbose && ps.verbosity > 2) {
      std::cout << "\tTotal number of cuts are " << total_num_cuts << " and " << total_num_legal_cuts << " were legal\n";
      std::cout << "\tSelected cut is " << max_i_1 << " " << max_i_2 << "\n";
      print_cut_cost (i, i_child[0], i_child[1], max_i_1, max_i_2, min_cost); 
    }
    assert((max_i_1 >= 0) && (max_i_2 >= 0));

    insert_cut_to_carry_list(i, ic, i_child[0], i_child[1], max_i_1, max_i_2); 

    return 0;
  }

  bool check_xilinx_5lut_legality ( uint32_t i1, uint32_t i2, uint32_t cut1, uint32_t cut2 ) {
    bool match = false;

    uint32_t nshared = 0;
    for (auto child_leaf1: cuts.cuts(i1)[cut1]) {
      for (auto child_leaf2: cuts.cuts(i2)[cut2]) {
        if (child_leaf1 == child_leaf2) {
          nshared++;
        } 
      }
    }
    if (cuts.cuts(i1)[cut1].size() + cuts.cuts(i2)[cut2].size() - nshared <= 5)
      return true;

    return false;
  }

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

      if (ps.verbose && ps.verbosity > 2) {
        std::cout << "Carry LUT mapping iteration " << i << ": " << n_carryin << " " << n_first << " " << n_second << "\n";
      }

      compute_carry_LUT_mapping(n_first, n_second, n_carryin);
    }
  }

  // This function takes 2 index with its carry-in index and figures out 
  // if it cuts can be mapped to LUTs before carry node and the cost
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
        double cost = cost_of_cuts(i1, i2, i1_child[0], i1_child[1],
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
          //print_cuts (cut_i1_1, cut_i1_2, cut_i2_1, cut_i2_2, i1_child, i2_child, cost); 
          print_cut_cost (i1, i1_child[0], i1_child[1], cut_i1_1, cut_i1_2, cost); 
          print_cut_cost (i2, i2_child[0], i2_child[1], cut_i2_1, cut_i2_2, cost); 
        }
          }
        }
      }
    }

    if (ps.verbose && ps.verbosity > 3) {
      std::cout << "\tTotal number of cuts are " << total_num_cuts << " and " << total_num_legal_cuts << " were legal\n";
      std::cout << "\tSelected cut is " << max_i1_1 << " " << max_i1_2 << " " << max_i2_1 << " " <<  max_i2_2 << "\n";
      print_cut_cost (i1, i1_child[0], i1_child[1], max_i1_1, max_i1_2, max_cost); 
      print_cut_cost (i2, i2_child[0], i2_child[1], max_i2_1, max_i2_2, max_cost); 
      //print_cuts (max_i1_1, max_i1_2, max_i2_1, max_i2_2, i1_child, i2_child, max_cost); 
          //print_cuts (cut_i1_1, cut_i1_2, i1_child[0], i1_child[1], cost); 
          //print_cuts (cut_i2_1, cut_i2_2, i2_child[0], i2_child[1], cost); 
    }
    assert((max_i1_1 >= 0) && (max_i1_2 >= 0) && (max_i2_1 >= 0) && (max_i2_2 >= 0));

    insert_cut_to_carry_list(i1, ic, i1_child[0], i1_child[1], max_i1_1, max_i1_2); 
    insert_cut_to_carry_list(i2, i1, i2_child[0], i2_child[1], max_i2_1, max_i2_2); 

    return 0;
  } 

  void update_delay () {

    for ( auto const& n : top_order ) {
      const auto index = ntk.node_to_index( n );

      if (ntk.is_pi(n) || ntk.is_constant(n)) {
        continue;
      } else if (is_a_carry_node(n)) {
        uint32_t updated_delay = 0;
        uint32_t updated_delay_index = 0;
        for (auto leaf: carry_cut_list[index]) {
          uint32_t leaf_delay = delays[leaf];
        
          // leafs from carry node get lut + adder delay
          if (carry_driver_nodes[index] == leaf) leaf_delay += CARRY_DELAY;
          else  leaf_delay += LUT_ADDER_DELAY; 

          if (leaf_delay > updated_delay) {
            updated_delay = leaf_delay;
            updated_delay_index = leaf;
          }
        }
        if (updated_delay > delays[index]) std::cout << "Carry delay has changed for node "
          << index << " by " << updated_delay_index << ": " << delays[index] << " -> " << updated_delay << "\n";
        if (updated_delay < delays[index]) std::cout << "Carry delay has changed for node "
          << index << ": " << delays[index] << " <- " << updated_delay << "\n";

        delays[index] = updated_delay;
      } else { 
        uint32_t new_delay = cut_flow(cuts.cuts(index).best()).second;
        if (new_delay > delays[index]) std::cout << "LUT delay has changed for node "
          << index << ": " << delays[index] << " -> " << new_delay << "\n";
        if (new_delay < delays[index]) std::cout << "LUT delay has changed for node "
          << index << ": " << delays[index] << " <- " << new_delay << "\n";
        delays[index] = new_delay;
      }
    } 
    delay = 0;
    ntk.foreach_po( [this]( auto s ) {
      const auto index = ntk.node_to_index( ntk.get_node( s ) );
      delay = std::max( delay, delays[index] );
    });


    print_all_best_cut();
  }

  // This function updates the depth information after carry mapping
  // Used for the next iteration of carry mapping
  void set_carry_mapping_refs()
  {
 
    /* compute current delay */
    delay = 0;
    ntk.foreach_po( [this]( auto s ) {
      const auto index = ntk.node_to_index( ntk.get_node( s ) );
      delay = std::max( delay, delays[index] );
      map_refs[index]++;
    });

    /*for (auto n: path_for_carry_chain) {
      map_refs[ntk.node_to_index(n)] = 0; // NOT NECESSARY?
      for (auto const cut_index: carry_cut_list[ntk.node_to_index(n)]) {
        if (!is_a_carry_node(ntk.index_to_node(cut_index)) && \
           !ntk.is_pi(ntk.index_to_node(cut_index))) {
          map_refs[cut_index]++;
        }
      }
    }*/

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
      
      if (is_a_carry_node( *it )) {
        map_refs[index]++;
        for ( auto leaf : carry_cut_list[index] )
        {
          //if (!is_a_carry_node(leaf))
          map_refs[leaf]++;
        }
  
      } else {
        for ( auto leaf : cuts.cuts( index )[0] )
        {
          map_refs[leaf]++;
        }
      }

      area++;
    }
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

  // This function gets the 2 children node that are driven from the LUT
  // Used for mapping the LUTs before carry
  void get_children_node (node<Ntk>child_index[2], uint32_t index, uint32_t carryin) {

    uint32_t curr_LUT = 0;
    auto const& curr_node = ntk.index_to_node(index);

    assert (ntk.fanin_size(curr_node) == 3);

    for (uint32_t i_fanin = 0; i_fanin < ntk.fanin_size(curr_node); i_fanin++) { 
      node<Ntk> leaf_node = ntk.get_children (curr_node, i_fanin);  
      
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
  uint64_t get_bit (uint64_t val, uint64_t pos) {
    return (val >> pos)&0x1;
  }
  uint64_t get_not_bit (uint64_t val, uint64_t pos) {
    if (((val >> pos)&0x1) == 1) return 0;
    return 0x1;
  }

  kitty::dynamic_truth_table compute_carry_function (uint32_t a, uint32_t b, uint32_t cindex_1, uint32_t cindex_2, bool child_complement[3], uint32_t unique_cuts[5], uint32_t n_total, 
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

    auto cut1_function_bu = cuts.truth_table( cut1 ); 
    auto cut2_function_bu = cuts.truth_table( cut2 ); 
    
    auto cut1_function = flip_carry_user_tt (cindex_1, cindex_1, a );
    auto cut2_function = flip_carry_user_tt (cindex_2, cindex_2, b );

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
    std::cout << "\tcut1old: ";
    kitty::print_hex(cut1_function_bu);
    std::cout << "\n\tcut2old: ";
    kitty::print_hex(cut2_function_bu);
    std::cout << "\n\tcut1: ";
    kitty::print_hex(cut1_function);
    std::cout << "\n\tcut2: ";
    kitty::print_hex(cut2_function);
    std::cout << "\n\tcombined: ";
    kitty::print_hex(function);
    std::cout << "\n";
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

  void determine_child_complement (bool* child_complement, node<Ntk> n, \
      uint32_t cindex_c, uint32_t cindex_1, uint32_t cindex_2) {
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
    } );
  }

  // This function takes the index and its children index with seleted best cut for
  // LUT before carry node.
  //  convention is that carry_cut_list contains 
  //  the inputs to the carry LUT (so push_back)*
  //  *maybe change to keep cut numbers
  // can't just keep cut numbers since we need to run cut enumeration...
  void insert_cut_to_carry_list(uint32_t index, uint32_t cindex_c,
      uint32_t cindex_1, uint32_t cindex_2, int32_t i, int32_t j) {

    // If index is 0, there is no cut for it
    // This is due to double LUT mapping that happens in ALMs
    // Only happens in the odd case for the last node
    if (index == 0) return;

    auto n = ntk.index_to_node(index);

    // Clear any existing cuts
    // Happens when we rerun after inverter removal
    carry_cut_list[index].clear();

    cut_t const& cut_1 = cuts.cuts(cindex_1)[i];
    cut_t const& cut_2 = cuts.cuts(cindex_2)[j];

    // Update delay data
    // cut_delay returns the delay of LUT
    // add LUT_DELAY to find the delay post carry node
    // NOT using LUT_ADDER_DELAY here (LUT_DELAY + LUT_DELAY)
    uint32_t cut_1_delay = cut_delay( cut_1 ) + LUT_ADDER_DELAY;
    uint32_t cut_2_delay = cut_delay( cut_2 ) + LUT_ADDER_DELAY;
    uint32_t carry_delay = delays[cindex_c] + CARRY_DELAY;
    uint32_t cut_delay = std::max( cut_1_delay, cut_2_delay);
    uint32_t new_delay = std::max( cut_delay, carry_delay );
    std::cout << index << ": " << cut_1 << "(" << cut_1_delay << "),"
                               << cut_2 << "(" << cut_2_delay << "),"
                               << cindex_c << "(" << carry_delay << ")\n";
    std::cout << "Delay changed for node " << index << ": " << delays[index] << " -> " << new_delay << "\n";  
    delays[index] = new_delay;

    // Add ONE of the LUTs to exit from ALM directly
    if (ps.lut_sharing) {
      if (cut_1.size() > 1 && mapped_to_5LUT[cindex_1] == false && (cut_2.size() <= 1 || ntk.fanout_size(cindex_1) > ntk.fanout_size(cindex_2)) ){
        mapped_to_5LUTa[index] = cindex_1;
        mapped_to_5LUT [ cindex_1 ] = true;
        std::cout << "Adding 1 " << cindex_1 << " as the LUT to exit\n";
      } else if (cut_2.size() > 1 && mapped_to_5LUT[cindex_2] == false  && (cut_1.size() <= 1 || ntk.fanout_size(cindex_2) > ntk.fanout_size(cindex_1)) ){
        mapped_to_5LUTb[index] = cindex_2;
        mapped_to_5LUT [ cindex_2 ] = true;
        std::cout << "Adding 2 " << cindex_2 << " as the LUT to exit\n";
      }
    }

    // Set LUT function
    // cut_1 and cut_2 with MIG carry node should be
    // turned into truth table 
    uint32_t index_child_cuts[2][5] = {0}; 
    uint32_t unique_cuts[5] = {0};
    uint32_t n_total = 0;   
    
    for (auto leaf: cut_1) {
      unique_cuts[n_total] = leaf;
      n_total++;
    }
    for (auto leaf: cut_2) {
      bool inserted = false;
      for (uint32_t j = 0; j < 5; j++) {
        if (unique_cuts[j] == leaf) inserted = true;
      }
      if (!inserted) {
        unique_cuts[n_total] = leaf;
        n_total++;
      }
    }

    // carry in node MUST be the last in the list for BLIF generation
    for (uint32_t i = 0; i < n_total; i++) {
      carry_cut_list[index].push_back(unique_cuts[i]);
    }
    carry_cut_list[index].push_back(cindex_c);
    carry_cut_index_list[index].push_back(i);
    carry_cut_index_list[index].push_back(j);

  }

  // Find the delay of the specific cut
  // It should look for the slowest (max) of all nodes
  uint32_t cut_delay ( cut_t const& cut ){
    uint32_t time{1u};

    for ( auto leaf : cut )
    {
      time = std::max( time, delays[leaf] );
    }

    return time;
  }

  // Lowest cost will get chosen!
  double cost_of_5LUT_cuts(uint32_t index, uint32_t index_1, uint32_t index_2, int32_t i1, uint32_t i2) {

    double cost = 0;
      
    if (ps.cost == 0) {
      cost += cuts.cuts(index_1)[i1].size() + cuts.cuts(index_2)[i2].size();

    } else if (ps.cost == 1) { // Best Delay
  
      cost += (1.0/(double)cut_delay (cuts.cuts(index_1)[i1]));
      cost += (1.0/(double)cut_delay (cuts.cuts(index_2)[i2]));

    } else if (ps.cost == 2) { // Most Number of Shared Inputs
  
      cost += 0;
    
    } else if (ps.cost == 3) { // Delay with shared input

      cost += (1.0/(double)cut_delay (cuts.cuts(index_1)[i1]));
      cost += (1.0/(double)cut_delay (cuts.cuts(index_2)[i2]));

    } else if (ps.cost == 4) {

     
      cost += (double)( cut_delay(cuts.cuts(index_1)[i1]) + LUT_DELAY) / delay;
      cost += (double)( cut_delay(cuts.cuts(index_2)[i2]) + LUT_DELAY) / delay;
      //cost += (double)((cuts.cuts(index_1)[i1].size() + cuts.cuts(index_2)[i2].size()))/10;
      
    } else {
      assert(0);
    }
    return cost;

  }

  double cost_of_cuts(uint32_t index1, uint32_t index2,
      uint32_t index1_1, uint32_t index1_2,
      uint32_t index2_1, uint32_t index2_2, uint32_t i1, uint32_t i2, uint32_t j1, uint32_t j2,
      uint32_t nshared1, uint32_t nshared2, uint32_t ntotal1, uint32_t ntotal2, uint32_t nshared) {

    double cost = 0;
  
    int cost_type = ps.cost;
    
    if (cost_type == 0) {
      cost += cuts.cuts(index1_1)[i1].size() + cuts.cuts(index1_2)[i2].size();
      cost += cuts.cuts(index2_1)[j1].size() + cuts.cuts(index2_2)[j2].size();

    } else if (cost_type == 1) { // Best Delay
  
      cost += (1.0/(double)cut_delay (cuts.cuts(index1_1)[i1]));
      cost += (1.0/(double)cut_delay (cuts.cuts(index1_2)[i2]));
      cost += (1.0/(double)cut_delay (cuts.cuts(index2_1)[j1]));
      cost += (1.0/(double)cut_delay (cuts.cuts(index2_2)[j2]));

    } else if (cost_type == 2) { // Most Number of Shared Inputs
  
      cost += 2*nshared; 
      cost += nshared1;
      cost += nshared2;
    
    } else if (cost_type == 3) { // Delay with shared input

      cost += (1.0/(double)cut_delay (cuts.cuts(index1_1)[i1]));
      cost += (1.0/(double)cut_delay (cuts.cuts(index1_2)[i2]));
      cost += (1.0/(double)cut_delay (cuts.cuts(index2_1)[j1]));
      cost += (1.0/(double)cut_delay (cuts.cuts(index2_2)[j2]));
      cost += 0.5*(nshared + nshared + nshared2); 

    } else if (cost_type == 4) {
      cost += (1.0/(double)(cut_delay (cuts.cuts(index1_1)[i1] )+ LUT_DELAY));
      cost += (1.0/(double)(cut_delay (cuts.cuts(index1_2)[i2] )+ LUT_DELAY));
      cost += (1.0/(double)(cut_delay (cuts.cuts(index2_1)[j1] )+ LUT_DELAY));
      cost += (1.0/(double)(cut_delay (cuts.cuts(index2_2)[j2] )+ LUT_DELAY));
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
    for ( auto const& n : top_order )
    {
      if ( ntk.is_constant( n ) || ntk.is_pi( n ) || is_a_carry_node( n ) )
        continue;
      compute_best_cut<ELA>( ntk.node_to_index( n ), delay );
    }
    set_mapping_refs<ELA>();
  }

  template<bool ELA>
  void set_mapping_refs()
  {

    const auto coef = 1.0f / ( 1.0f + ( iteration + 1 ) * ( iteration + 1 ) );

    /* compute current delay and update mapping refs */
    delay = 0;
    ntk.foreach_po( [this]( auto s ) {
      const auto index = ntk.node_to_index( ntk.get_node( s ) );
      delay = std::max( delay, delays[index] );

      if constexpr ( !ELA ) {
        map_refs[index]++;
      }
    } );


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
      flow += flows[leaf];
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

      map_refs[leaf]--;
      // This should never be negative (unsigned integer)
      assert(map_refs[leaf] <1000000);
      if ( map_refs[leaf] == 0 )
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
        time = cut_flow( *cut ).second;
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
      if (best_time > delays[index]) std::cout << "ELA changed " << index 
        << ": from " << delays[index] << "->" << best_time << "\n";
      if (best_time < delays[index]) std::cout << "ELA changed " << index 
        << ": from " << delays[index] << "<-" << best_time << "\n";
    }
    delays[index] = best_time;
    flows[index] = best_flow / flow_refs[index];
    
    if ( best_cut != 0 )
    {
      cuts.cuts( index ).update_best( best_cut );
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
        auto const n_carry = path[i-1];

        if ( ntk.is_constant( n ) || ntk.is_pi( n ) || ntk.is_cell_root( n ) )  
          continue;

        if ( map_refs[n] == 0 ) assert(0);
        
        const auto index = ntk.node_to_index( n );

        if (ps.verbose && ps.verbosity > 3) std::cout << "Node* " << n;

        // Add node as carry
        std::vector<node<Ntk>> nodes;
        ntk.add_to_carry_mapping ( n, ntk.node_to_index(n_carry) ); 
        
        // Add LUT mapped before carry  
        uint32_t index_c1 = 0;
        uint32_t index_c2 = 0;
        if (!is_a_carry_node(mapped_to_5LUTa[index]) && mapped_to_5LUTa[index] > 0) 
          index_c1 = mapped_to_5LUTa[index];
        if (!is_a_carry_node(mapped_to_5LUTb[index]) && mapped_to_5LUTb[index] > 0)
          index_c2 = mapped_to_5LUTb[index];

        ntk.add_to_carry_LUT_mapping( n, index_c1, index_c2, \
            mapped_to_5LUT_complemented[index_c1], mapped_to_5LUT_complemented[index_c2] );
        if (ps.verbose && ps.verbosity > 3) std::cout << "(" << index_c1 << "/" << index_c2 << ")" << ":";

        for ( auto c: carry_cut_list[ntk.node_to_index(n)]) {
          nodes.push_back(c);
          if (ps.verbose && ps.verbosity > 3) std::cout << c << " ";
        }
        if (ps.verbose && ps.verbosity > 3) std::cout << "\n";
        if (ps.verbose && ps.verbosity > 3) std::cout << "\t";
        if (ps.verbose && ps.verbosity > 3) kitty::print_hex(ntk.cell_function(n));
        if (ps.verbose && ps.verbosity > 3) std::cout << "\n";
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
      if ( ntk.is_constant( n ) || ntk.is_pi( n ) || ntk.is_cell_root( n ) || is_a_carry_node( n ) ) 
        continue;

      const auto index = ntk.node_to_index( n );

      if ( map_refs[index] == 0 ) 
        continue;
      if (ps.verbose && ps.verbosity > 3) std::cout << "Node " << n << ": ";
      std::vector<node<Ntk>> nodes;
      for ( auto const& l : cuts.cuts( index ).best() ) {
        nodes.push_back( ntk.index_to_node( l ) );
        if (ps.verbose && ps.verbosity > 3) std::cout << l << " ";
      }
      if (ps.verbose && ps.verbosity > 3) std::cout << "\n";
      ntk.add_to_mapping( n, nodes.begin(), nodes.end() );

      // This is translated into LUT functionality
      if constexpr ( StoreFunction )
      {
        if ( !is_a_carry_node(n) ) {
          // NEED TO CHANGE THIS
          auto tt_updated = flip_carry_user_tt ( index, index, 0 );

          ntk.set_cell_function( n, tt_updated );

          if (ps.verbose && ps.verbosity > 3) std::cout << "\t";
          if (ps.verbose && ps.verbosity > 3) kitty::print_hex(ntk.cell_function(n));
          if (ps.verbose && ps.verbosity > 3) std::cout << "\n";
        }
      }
    }
  }

  kitty::dynamic_truth_table flip_carry_user_tt ( uint32_t index, uint32_t curr_index, uint32_t cut_i ) {
 
    kitty::dynamic_truth_table tt_new (cuts.cuts(index)[cut_i].size());;

    if (ntk.is_constant(curr_index)) return tt_new;

    uint32_t place = 0;
    bool here = false;
    bool flip = false;
    for (auto leaf: cuts.cuts(index)[cut_i]) {
      if (leaf == curr_index ) {
        if (mapped_to_5LUT_complemented[leaf] == true) {
          std::cout << "\tflip " << leaf << "\n";
          flip = true;
        }
        here = true;
        break;
      }
      place++;
    } 
 
    if (here) {
      for (uint64_t k = 0; k < pow(2,cuts.cuts(index)[cut_i].size()); k++) {
        if (flip == true)
          tt_new._bits[0] += uint64_t(get_not_bit(k, place) << k);
        else 
          tt_new._bits[0] += uint64_t(get_bit(k, place) << k);
      }
      return tt_new;
    }

    auto n = ntk.index_to_node(curr_index);

    uint32_t j = 0; 
    std::vector<kitty::dynamic_truth_table> tt_child(3);
    ntk.foreach_fanin( n, [&]( auto const& f ) {
      auto n_child = ntk.get_node(f);
      tt_child[j] = flip_carry_user_tt(index, n_child, cut_i);
      j++;  
    });
  
    tt_new = ntk.compute(n, tt_child.begin(), tt_child.end());
    return tt_new;
    
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

  void print_inverter (void) {
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

  void print_cut_cost (uint32_t i, uint32_t i1, uint32_t i2, uint32_t cut1, uint32_t cut2, double cost) {
    std::cout << "\tCut " << cut1 << " " << cut2 <<": ";
    std::cout << cost << "\n";
    std::cout << "\t\t" << i1 << ":" << cuts.cuts(i1)[cut1] << " " << cut_delay(cuts.cuts(i1)[cut1]) + LUT_DELAY << "\n";
    std::cout << "\t\t" << i2 << ":" << cuts.cuts(i2)[cut2] << " " << cut_delay(cuts.cuts(i2)[cut2]) + LUT_DELAY << "\n";
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
        std::cout << leaf;
        if (delays[leaf] > worst_delay) {
          critical_index = leaf;
          worst_delay = delays[leaf];
        }
        std::cout << "(" << delays[leaf] << ") ";
      }
    } else {
      for (auto leaf: cuts.cuts(index)[0]) {
        std::cout << leaf;
        if (delays[leaf] > worst_delay) {
          critical_index = leaf;
          worst_delay = delays[leaf];
        }
        std::cout << "(" << delays[leaf] << ") ";
      }
    }
    std::cout << " -- " << critical_index <<"\n";
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
    std::cout << "total lut count and carry count: " << lut_count << " " << carry_count << "\n";
  }

  // Print best cuts of a node
  void print_all_best_cut() {

    for ( auto i = 0u; i < ntk.size(); ++i ) {

      uint32_t max_delay = 0;
      uint32_t max_index = 0;
      if (is_a_carry_node(i)) {
        std::cout << "Node* " << i << "(" << delays[i] << "): { ";
        for ( auto leaf: carry_cut_list[i] ){
          std::cout << leaf << " ";
          if ( delays[leaf] > max_delay ) {
            max_delay = delays[leaf];
            max_index = leaf;
          }
        }
        std::cout << "} -> " << max_index <<"\n";
      } else {
        std::cout << "Node " << i << "(" << delays[i] << "): " << cuts.cuts(i).best();
        for ( auto leaf: cuts.cuts(i).best() ){
          if ( delays[leaf] > max_delay ) {
            max_delay = delays[leaf];
            max_index = leaf;
          }
        }
        std::cout << " -> " << max_index <<"\n";
      }
    }
  }


  void print_state()
  {

    ntk.foreach_node( [&]( auto n, auto ) {
      if ( ntk.fanout_size(n) == 0 ) 
        std::cout << "Node " << n << " has 0 output\n";
    });
    ntk.foreach_po( [&]( auto const& s ) {
      std::cout << "PO " << ntk.get_node(s) << ":" << delays[ntk.get_node(s)]<<"\n";
    });
    if (ps.verbosity > 4 ) {
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
  //bool ela{false};       /* compute exact area */

  std::vector<node<Ntk>> top_order;
  std::vector<float> flow_refs;
  std::vector<uint32_t> map_refs;
  std::vector<float> flows;
  std::vector<uint32_t> delays;
  std::vector<std::vector<uint32_t>> carry_cut_list;
  std::vector<std::vector<uint32_t>> carry_cut_index_list;
  std::vector<node<Ntk>> carry_nodes;
  std::vector<node<Ntk>> carry_driver_nodes;
  std::vector<std::vector<node<Ntk>>> carry_paths;
  std::vector<node<Ntk>> mapped_to_5LUTa;
  std::vector<node<Ntk>> mapped_to_5LUTb;
  std::vector<bool> mapped_to_5LUT;
  std::vector<bool> mapped_to_5LUT_complemented;

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
