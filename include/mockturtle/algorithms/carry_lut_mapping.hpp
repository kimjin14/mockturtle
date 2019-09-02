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

#define CARRY_DELAY 1
#define LUT_DELAY 3

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
  bool verbose{false};

  /*! \brief Map to carry. */
  bool carry_mapping{true};

  /*! \brief Number of rounds for carry chain synthesis. */
  uint32_t rounds_carry{3u};  
  
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
        carry_cut_index_list ( ntk.size() ),
        carry_paths(ps.rounds_carry),
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
 
    if (ps.carry_mapping) {

      std::vector<node<Ntk>> path_for_carry;

      // Initial mapping to LUTs only.
      set_mapping_refs<false>();
      compute_mapping<false>();

      delay_lut = delay;

      for (uint32_t i = 0; i < ps.rounds_carry; i++) { 

        //std::cout << "Carry Mapping Iteration " << i << " targetting delay of " << delay_lut << "\n";

        /////////////////////////////////////////////////
        // Find path
        // find_critical_paths(path_for_carry);
        find_critical_LUT_chain (path_for_carry);
        //print_path (path_for_carry);
        carry_paths[i] = path_for_carry;
        for (uint32_t j = 1; j < path_for_carry.size(); j++) {
          carry_driver_nodes[path_for_carry[j]] = path_for_carry[j-1];
        }
        if (path_for_carry.size() == 0) break;

        ////////////////////////////////////////////////
        // Compute carry LUT mappgin
        //map_carry_lut_with_node(path_for_carry);
        compute_carry_mapping<false>(path_for_carry);
        set_mapping_refs<false>();
        compute_mapping<false>();

        ////////////////////////////////////////////////
        // Clear path 
        path_for_carry.clear();
      }
      
      // Must update cut truth tables since inverters were moved
      // Currently rerunning cut_enumeration but TODO: need to just update truth table
      // ////////////////////////////////////////////////
      // Remove inverters in its path
      remove_inverter();
      cuts = cut_enumeration<Ntk, StoreFunction, CutData>( ntk, ps.cut_enumeration_ps ); 
      for (auto update_carry_path: carry_paths)
        compute_carry_mapping<true>(update_carry_path);
      //print_carry_paths();

      check_inverter();
      print_state();
    }    

    set_mapping_refs<false>();

    while ( iteration < ps.rounds + (ps.rounds_carry*2 + 2) )
    {
      compute_mapping<false>();
    }

    while ( iteration < ps.rounds + ps.rounds_ela + (ps.rounds_carry*2 + 2) )
    {
      compute_mapping<true>();
    }
    derive_mapping();
    
  }

private:
  uint32_t cut_area( cut_t const& cut ) const
  {
    return static_cast<uint32_t>( cut->data.cost );
  }

  bool add_find_longest_path_in_mapping (std::vector<node<Ntk>>& path_for_carry, 
      uint32_t source_index, uint32_t dest_index) {

    auto const& n = ntk.index_to_node(source_index);
    
    ntk.set_visited(n, 1);

    if (source_index == dest_index) {
      return true;
    }
 
    bool found = false; 
    for (uint32_t i = 0; i < ntk.fanin_size(n); i++) {
      auto const leaf = ntk.get_children(n,i);
      auto const leaf_index = ntk.node_to_index(leaf);
      if (!is_a_carry_node(leaf) && !ntk.is_constant(leaf) && ntk.visited(leaf) == 0)
        found = add_find_longest_path_in_mapping (path_for_carry, leaf_index, dest_index);
      if (found) {
        path_for_carry.push_back(leaf);
        carry_nodes[leaf] = 1;
        //std::cout << "adding " << leaf << "\n";
        return true;
      }
    }
    return found;
  }

  // Keep going deeper to its fanin until depth is found
  bool find_deepest_LUT (std::vector<node<Ntk>>& path_for_carry, uint32_t index) { 
    //std::cout << "Index " << index << "(" << delays[index] << "):\n " ; 


    if (delays[index] == 0) {
      return true;
    } 

    bool deepest = false;

    uint32_t max_leaf = index;
    uint32_t max_delay = 0;

    if (is_a_carry_node(ntk.index_to_node(index))) {
      for ( auto leaf : carry_cut_list[index] ) {
        if (!is_a_carry_node(leaf)) { 
          //std::cout << "\tcarry " << leaf << "(" << delays[leaf] << ")\n"; 
          if (max_delay < delays[leaf]) {
            max_leaf = leaf;
            max_delay = delays[leaf];
          }
        }
      }
    } else {
      for ( auto leaf : cuts.cuts( index )[0] ){
        if (!is_a_carry_node(leaf)) { 
          //std::cout << "\t\treg " << leaf << "(" << delays[leaf] << ")\n"; 
          if (max_delay <= delays[leaf]) {
            max_leaf = leaf;
            max_delay = delays[leaf];
          }
        }
      }
    }

    if (max_leaf != index) {
      deepest = find_deepest_LUT(path_for_carry, max_leaf);
      if (deepest)  {
        ntk.clear_visited();
        //std::cout << index << " to " << max_leaf << "\n";
        if(add_find_longest_path_in_mapping(path_for_carry, index, max_leaf)) {
          //std::cout << "add the path\n";
        } else { //cannot put this in carry
          //std::cout << "what should happen here?\n";
          return false;
        } 
  
      }
    }
    return deepest;
  }
  
  // For each output, recursively go to the next node
  //  until the deepest LUT depth is found
  void find_critical_LUT_chain(std::vector<node<Ntk>>& path_for_carry) {

    // Search through the nodes in reverse order 
    for ( auto it = top_order.rbegin(); it != top_order.rend(); ++it ) {
      const auto index = ntk.node_to_index( *it );
      if (delays[index] == delay_lut && delay_lut!= 0 && !ntk.is_pi(*it)) {
        //std::cout << "\nFound target index " << index << "\n";
        if (find_deepest_LUT(path_for_carry,index)) {
          path_for_carry.push_back(index);
          carry_nodes[index] = 1;
        }
        break;
      }
    }
  }

  // These functions find a path that is not already on the carry chain
  bool get_path(std::vector<node<Ntk>>& path_for_carry, 
      node<Ntk> n, uint32_t depth, uint32_t curr_depth) {

    ntk.set_visited(n, 1);

    if (ntk.is_pi(n)) {
      if (curr_depth == depth) {
        carry_nodes[n] = 1;
        path_for_carry.push_back(n);
        return true;
      } else return false;
    }

    bool longest = false;

 
    for (uint32_t i = 0; i < ntk.fanin_size(n); i++) {
      auto nchild = ntk.get_children(n,i);
      if (is_a_carry_node(nchild) || is_in_carry_lut(nchild) || ntk.visited(nchild)!=0 ) continue;
      longest = get_path(path_for_carry, nchild, depth, curr_depth+1);
      if (longest) {
        carry_nodes[n] = 1;
        path_for_carry.push_back(n);
        break;
      }
    }
    return longest;
  }

  void find_critical_paths(std::vector<node<Ntk>>& path_for_carry) {

    auto depth = depth_view<Ntk>(ntk).depth();

    for (uint32_t i = 0; i < ntk.num_pos(); i++) {
      auto n = ntk.get_po(i);
      if (is_a_carry_node(n) || is_in_carry_lut(n)) continue;

      ntk.clear_visited();
      if (get_path(path_for_carry,n,depth,0)) {
        break;
      }
    }
  }
  // If the node is a carry, check its carry child for inversion
  void check_inverter (void) {

    for (auto carry_path: carry_paths) {
      for (uint32_t carry_i = 1; carry_i < carry_path.size(); carry_i++) {

        auto n = carry_path[carry_i];
        auto n_carry = carry_path[carry_i-1];

        //std::cout << "node " << n << ": ";
        for (uint32_t i = 0; i < ntk.fanin_size(n); i++) {

          node<Ntk> child_node = ntk.get_children(n,i);  
          //std::cout << child_node << "(" << ntk.is_complemented_children(n,i) << ")";
          if(child_node == n_carry) {
            //std::cout << "*";
            if (ntk.is_complemented_children(n,i)) 
              assert(0);
          }
          //std::cout << " ";
        } 
        //std::cout << "\n";
      }
    } 
  }

  void remove_inverter ( ) {
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
            //std::cout << "\t\tflipping child " << child_node << "\n";
            ntk.flip_children(carry_node, i);
          }
          
          ntk.foreach_node( [&]( auto n, auto ) {
            for (uint32_t i = 0; i < ntk.fanin_size(n); i++) {
              if (ntk.get_children(n,i) == carry_node) {
                //std::cout << "\t\tflipping child " << ntk.get_children(n,i) << \
                    " of node " << n << "\n";
                //std::cout << "\t\t\tfrom " << ntk.is_complemented_children(n,i); 
                ntk.flip_children(n, i);
                //std::cout << " to " << ntk.is_complemented_children(n,i) << "\n"; 
              }
            }
          });
          ntk.foreach_po( [&]( auto const& s ) {
            if (ntk.get_node(s) == carry_node) {
              //std::cout << "\t\tflipping output " << ntk.get_node(s) << "\n";
              ntk.flip_complement_output(s);
            }
          });

        }
      }
    }
  }

  void print_inverter_from_path (void) {

    //for (auto n: path_for_carry) {
    ntk.foreach_node( [this]( auto n, auto ) {
      std::cout << n << ":"; 
      for (uint32_t i = 0; i < ntk.fanin_size(n); i++) {
        node<Ntk> child_node = ntk.get_children(n,i);  
        std::cout << "  " << child_node << "(" << ntk.is_complemented_children(n,i) << ") ";
      }
      std::cout << "\n";
    });
  }

  void init_nodes()
  {
    ntk.foreach_node( [this]( auto n, auto ) {
      const auto index = ntk.node_to_index( n );
      map_refs[index] = 0;

      if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
      {
        /* all terminals have flow 1.0 */
        flow_refs[index] = 1.0f;
      }
      else
      {
        flow_refs[index] = static_cast<float>( ntk.fanout_size( n ) );
      }

      flows[index] = cuts.cuts( index )[0]->data.flow;
      delays[index] = cuts.cuts( index )[0]->data.delay;
    } );
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

  // Decide which nodes can fit into LUT
  // Set mapping ref for those nodes
  void map_carry_lut_with_node(std::vector<node<Ntk>> path_for_carry)
  {
    uint32_t starting_index = 0;
    if (ntk.is_pi(path_for_carry[0]) || is_a_carry_node(path_for_carry[0]))
      starting_index = 1;

    // Map nodes following the critical path in pairs (due to ALM)
    for (uint32_t i = starting_index; i < path_for_carry.size(); i+=2) {
      auto& n_carryin = path_for_carry[i-1];
      auto& n_first = path_for_carry[i];
      auto& n_second = path_for_carry[i+1];

      if (i == path_for_carry.size()-1) // odd number of nodes
        check_child_node(n_first, 0, n_carryin);
      else
        check_child_node(n_first, n_second, n_carryin);
    }
  }

  // Decide which nodes can fit into LUT
  // Set mapping ref for those nodes
  template<bool SetLUT>
  void compute_carry_mapping (std::vector<node<Ntk>> path_for_carry)
  {

    uint32_t starting_index = 0;
    node<Ntk> n_carryin;
    node<Ntk> n_first;
    node<Ntk> n_second;

    if (!path_for_carry.empty() && (ntk.is_pi(path_for_carry[0]) || is_in_carry_lut(path_for_carry[0]))) {
      starting_index = 1;
    }

    // Map nodes following the critical path in pairs (due to ALM)
    for (uint32_t i = starting_index; i < path_for_carry.size(); i+=2) {

      if (i == 0) n_carryin = 0;
      else n_carryin = path_for_carry[i-1];

      n_first = path_for_carry[i];
      n_second = path_for_carry[i+1];

      //std::cout << "i " << i << ": " << n_carryin << " " << n_first << " " << n_second << "\n";

      if (i == path_for_carry.size()-1) {// odd number of nodes
        n_second = 0;
      }

      if (n_second != 0 && (is_in_carry_lut(n_first) && !is_in_carry_lut(n_second))) {
        //std::cout << n_first << "(" << is_in_carry_lut(n_first) << ") " << \
          n_second << "(" << is_in_carry_lut(n_second) << ")\n";  
        i++;
        n_first = path_for_carry[i];
        n_first = path_for_carry[i+1];
        assert(0);
      } else if (n_second != 0 && (!is_in_carry_lut(n_first) && is_in_carry_lut(n_second))) {
        //std::cout << n_first << "(" << is_in_carry_lut(n_first) << ") " << \
          n_second << "(" << is_in_carry_lut(n_second) << ")\n";  
        n_second = 0;
        i--;
        assert(0);
      }
      check_cut_carry<SetLUT>(n_first, n_second, n_carryin);
    }
  }

  bool is_a_carry_node(node<Ntk> n) {
    if (carry_nodes[n] != 0) return true;
    return false;
  }

  bool is_in_carry_lut(node<Ntk> n) {
    if (!carry_cut_list[ntk.node_to_index(n)].empty()) return true;
    return false;
  }

  template<bool ELA>
  void compute_mapping()
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
        for ( auto leaf : carry_cut_list[ntk.node_to_index(n)] )
        {
          //time = std::max( time, delays[leaf] );
          if (delays[leaf] > time) {
            time = delays[leaf];
            slowest_leaf = leaf;
          }
        }
        if (carry_driver_nodes[n] == slowest_leaf)
          delays[ntk.node_to_index(n)] = time + CARRY_DELAY;
        else 
          delays[ntk.node_to_index(n)] = time + LUT_DELAY;
        continue;
      }  
      compute_best_cut<ELA>( ntk.node_to_index( n ) );
    }
    set_mapping_refs<ELA>();
  }

  template<bool ELA>
  void set_mapping_refs()
  {

    //std::cout << "set_mapping_refs " << iteration << "\n";
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
            delay_lut = std::max( delay_lut, delays[cut_index]);
            map_refs[cut_index]++;
          }
        }
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

  // a b c e0
  // a b c f0
  // 5 unique inputs are possible, 3 shared
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

  bool check_5lut_legality (uint32_t input_array[2][5], \
      uint32_t unique_array[8], uint32_t* n_shared, uint32_t* n_total) {

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

    // This makes sure that we are dealing with MIG nodes
    //std::cout << curr_node << "(" << ntk.fanin_size(curr_node) << ")\n";
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
  void list_covered_nodes(node<Ntk> n, cut_t const& cut) {

    if (cut.size() == 1) {
      for (auto const &l : cut) {
        if (l == ntk.node_to_index(n)) {
          return;
        } else {
          assert(0);
        }
      }
    } 

    ntk.set_visited(n,1);
    
    carry_cut_list[ntk.node_to_index(n)].push_back(1);

    for (uint32_t i = 0; i < ntk.fanin_size(n); i++) {
      bool found = false;
      auto c = ntk.get_children(n,i);
      for ( auto const& l : cut ) {
        if (l == ntk.node_to_index(c)) {
          found = true;
        } 
      } 
      if (!found && !ntk.is_pi(c) && !ntk.is_constant(c) && !(ntk.visited(c)>0)) {
        carry_cut_list[ntk.node_to_index(c)].push_back(1);
        list_covered_nodes(c, cut);
      }
    }
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

    //std::cout << "\tUnique list(" << n_total << "):";
    //for (uint32_t i = 0; i < n_total; i++) {
      //std::cout << unique_cuts[i] << ",";
    //}
    //std::cout << "\n";
    
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

    //std::cout << "cut1_index:";
    //for (uint32_t i = 0; i < n_total; i++) {
    //  std::cout << cut1_index[i] << ",";
    //}
    //std::cout << "\n";
    //std::cout << "cut2_index:";
    //for (uint32_t i = 0; i < n_total; i++) {
    //  std::cout << cut2_index[i] << ",";
    //}
    //std::cout << "\n";
    
    // Create a truth table
    uint32_t lut_input_size = n_total + 1; 
    kitty::dynamic_truth_table function (lut_input_size);
    function._bits[0]= 0x0000000000000000;

    auto cut1_function = cuts.truth_table( cut1 ); 
    auto cut2_function = cuts.truth_table( cut2 ); 
    
    //std::cout << "\tCut 1 Function:"; 
    //kitty::print_hex(cut1_function);
    //std::cout << "\n";
    //std::cout << "\tCut 2 Function:"; 
    //kitty::print_hex(cut2_function);
    //std::cout << "\n";
    

    // Go through each combination of the input
    // i4 i3 i2 i1 i0
    // 0 0 0 0 0
    // 0 0 0 0 1
    //std::cout << "\tFunction combination\n";
    uint32_t n_combinations = pow(2,n_total+1);
  
    if (ps.carry_lut_combined) {
      for (uint32_t i = 0; i < n_combinations; i++) {
        /*std::cout << "\t";
        for (int32_t j = n_total; j >= 0; j--) {
          std::cout << get_bit(i, j) << " ";
        } 
        std::cout << "\n";
        */
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
        
        //std::cout << "\t" << cut1_val << " " << cut2_val << " " << carry_val << " -> " << lut_val << "\n";
        //  Set that value in the right place of the function
        function._bits[0] += uint64_t( lut_val << i);
        //kitty::print_hex(function);
        //std::cout << "\n";
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
        
        //std::cout << "\t" << cut1_val << " " << cut2_val << " " << carry_val << " -> " << lut_val << "\n";
        //  Set that value in the right place of the function
        function._bits[0] += uint64_t( lut_val << i);
        //kitty::print_hex(function);
        //std::cout << "\n";
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
        //std::cout << "child is " << ntk.get_node(f) << " and does not fit " << \
            cindex_1 << " " << cindex_2 << " " << cindex_c << "\n";
        assert(0 && "Children index do not match with fanin.");
      }
      //std::cout << ntk.get_node(f) << "(" << ntk.is_complemented(f) << "),";
    } );
    //std::cout << "\n";

    

    // Cannot do this due to immutable view
    //if (!ps.carry_lut_combined) {
      //std::cout << "node is " << n << " to ";
      //auto cnodes = ordered_children(cindex_1);      
      //auto csignal = ntk.create_maj(cnodes[0], cnodes[1], cnodes[2]);
      //ntk.substitute_node(cindex_1, csignal);
      //std::cout << ntk.get_node(csignal) << "\n";
    //}

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
        //std::cout << "\tCannot map any cuts, need to just put the node\n";
        //std::cout << "\tfrom mig store function 2 " << n << ": ";
        //kitty::print_hex(mig_function);
        //std::cout << "\n";
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
      //std::cout << "\tfrom mig store function 1 " << n << ": ";
      //kitty::print_hex(function);
      //std::cout << "\n";
    }

    for (uint32_t i = 0; i < n_total; i++) {
      carry_cut_list[index].push_back(unique_cuts[i]);
    }
    carry_cut_list[index].push_back(cindex_c);
  }

  // This function takes 2 index with its carry-in index and figurse out 
  // if it cuts can be mapped to LUTs before carry node
  // Check if inputs are legal
  // each index can only have 5 unique inputs
  // index1 can't have more than 5 unique inputs
  //    e0 c b a
  //    f0 c b a
  //    e1 d b a
  //    f1 d b a
  template <bool SetLUT>
  int check_cut_carry(uint32_t index1, uint32_t index2, uint32_t indexc) {

    int32_t cost = 0;
    int32_t max_cost = 4;
    int32_t max_i = -1;
    int32_t max_j = -1;
    int32_t max_k = -1;
    int32_t max_l = -1;

    //std::cout << "index " << index1 << " " << index2 << " with carry " << indexc << "\n";

    node<Ntk> index1_child[2] = {0};
    node<Ntk> index2_child[2] = {0};

    get_children_node (index1_child, index1, indexc); 
    if (index2 != 0) get_children_node (index2_child, index2, index1); 

    for (uint32_t cut_index1_1 = 0; cut_index1_1 < cuts.cuts(index1_child[0]).size(); cut_index1_1++) {
      for (uint32_t cut_index1_2 = 0; cut_index1_2 < cuts.cuts(index1_child[1]).size(); cut_index1_2++) {
        for (uint32_t cut_index2_1 = 0; cut_index2_1 < cuts.cuts(index2_child[0]).size(); cut_index2_1++) {
          for (uint32_t cut_index2_2 = 0; cut_index2_2 < cuts.cuts(index2_child[1]).size(); cut_index2_2++) {

        cost = 0;


        // Array holding the nodes to the halfs of ALM as separate 4-LUT
        uint32_t index1_child_cuts[2][5] = {0}; 
        uint32_t index2_child_cuts[2][5] = {0}; 

        uint32_t abce0f0[8] = {0};
        uint32_t abde1f1[8] = {0};
  
        uint32_t ntotal_index1 = 0;
        uint32_t ntotal_index2 = 0;
   
        // compare first 2 4-LUT and count matches and unique innputs 
        uint32_t nshared_index1 = 0;
        uint32_t nshared_index2 = 0;

        bool index1_status = false;
        bool index2_status = false;

        uint32_t cost_index1 = 0;
        uint32_t cost_index2 = 0;

        // insert_cuts_to_array checks whether cut index can fit into the 4-LUT
        // check_5lut_legality checks whether those cuts can map to 5-LUT
        if (insert_cuts_to_array (index1_child_cuts, index1, indexc, cut_index1_1, cut_index1_2)) {
          index1_status = check_5lut_legality (index1_child_cuts, abce0f0, &nshared_index1, &ntotal_index1);
        }
      
        // if index2 is 0, it's not used 
        if (index2 != 0) {
          if (insert_cuts_to_array (index2_child_cuts, index2, index1, cut_index2_1, cut_index2_2)) { 
            index2_status = check_5lut_legality (index2_child_cuts, abde1f1, &nshared_index2, &ntotal_index2);
          }
        } 

        if (index1_status && index2_status) { 
          // At this point, both 5 LUTs are legal
          // Now we need to check for 6 LUT input legality
          // Depending on the total number of inputs, 
          uint32_t nshared_matches = 0; 
          for (uint32_t i = 0; i < nshared_index1; i++) {
            assert (abce0f0[i] != 0);
            for (uint32_t j = 0; j < nshared_index2; j++) {
              if (abce0f0[i] == abde1f1[j]) nshared_matches++;
            }
          }
          uint32_t nmatches = 0; 
          for (uint32_t i = 0; i < 5; i++) {
            if (abce0f0[i] == 0) continue;
            for (uint32_t j = 0; j < 5; j++) {
              if (abce0f0[i] == abde1f1[j]) nmatches++;
            }
          }
          assert(nshared_index1 <= 3 && nshared_index2 <= 3);
          if (((ntotal_index1 + ntotal_index2 - nmatches) <= 8) &&
            (int(nshared_matches) >= int((nshared_index1 - 1) + (nshared_index2 - 1) - \
                (5 - ntotal_index1) - (5 - ntotal_index2)))) { 
          } else {
            if (cost_index1 > cost_index2) index2_status = false;
            else index1_status = false;
          }
        }

        if (index1_status) {
          cost_index1 = cuts.cuts(index1_child[0])[cut_index1_1].size() +
              cuts.cuts(index1_child[1])[cut_index1_2].size();
          cost += cost_index1;
        } else {
          cost += 2;
        }

        if (index2_status) {
          cost_index2 = cuts.cuts(index2_child[0])[cut_index2_1].size() +
              cuts.cuts(index2_child[1])[cut_index2_2].size();
          cost += cost_index2; 
        } else {
          cost += 2;
        }


        if (cost > max_cost) {
          //std::cout << cut_index1_1 << " " << cut_index1_2 << " " << cut_index2_1 << " " << cut_index2_2 << " -> " << cost << "\n";
          max_cost = cost;
          if (index1_status) {
            max_i = cut_index1_1;
            max_j = cut_index1_2;
          }
          if (index2_status) {
            max_k = cut_index2_1;
            max_l = cut_index2_2;
          }
        }

          }
        }
      }
    }

    //carry_cut_index_list[index1].push_back(max_i);
    //carry_cut_index_list[index1].push_back(max_j);
    //carry_cut_index_list[index2].push_back(max_k);
    //carry_cut_index_list[index2].push_back(max_l);

    insert_cut_to_carry_list<SetLUT>(index1, indexc, index1_child[0], index1_child[1], max_i, max_j); 
    insert_cut_to_carry_list<SetLUT>(index2, index1, index2_child[0], index2_child[1], max_k, max_l); 

    return 0;
  } 


  // Takes 2 indices and the carryin node value and determines whether its
  // children can fit into the LUTs preceding FA
  int check_child_node (uint32_t index1, uint32_t index2, uint32_t carryin) {

    std::cout << "mapping " << index1 << " " << index2 << "\n";

    // Array holding the nodes to the halfs of ALM as separate 4-LUT
    uint32_t index1_inputs[2][5] = {0}; 
    uint32_t index2_inputs[2][5] = {0}; 
    node<Ntk> index1_child[2] = {0};
    node<Ntk> index2_child[2] = {0};

    // Each node should have 3 children, where 1 is mapped to carry in
    // 2 inputs should use LUTs with above restrictions
    insert_child_to_array (index1_inputs, index1_child, index1, carryin);
    if (index2 != 0) {
      insert_child_to_array (index2_inputs, index2_child, index2, index1);
    }     
    // Check if inputs are legal
    // each index can only have 5 unique inputs
    // index1 can't have more than 5 unique inputs
    //    e0 c b a
    //    f0 c b a
    //    e1 d b a
    //    f1 d b a
    uint32_t abce0f0[8] = {0,0,0,0,0};
    uint32_t abde1f1[8] = {0,0,0,0,0};
  
    uint32_t n_unique_index1 = 0;
    uint32_t n_unique_index2 = 0;
   
    // compare first 2 4-LUT and count matches and unique innputs 
    uint32_t n_shared_index1 = 0;
    uint32_t n_shared_index2 = 0;
    n_unique_index1 = find_unique_set (&n_shared_index1, index1_inputs, abce0f0);
    n_unique_index2 = find_unique_set (&n_shared_index2, index2_inputs, abde1f1);

     // Add carry in to the carry cut list
    carry_cut_list[index1].push_back(carryin);
    if (index2!=0) carry_cut_list[index2].push_back(index1);

    // node1 and node 2 cannot fit (too many inputs)
    if (n_unique_index1 > 5 && n_unique_index2 > 5) {
      for (uint32_t i = 0; i < 2; i++) { 
        if (index1_child[i] != 0) {
          carry_cut_list[index1].push_back(index1_child[i]);
          //carry_cut_list[index1_child[i]].push_back(index1);
        }
        if (index2_child[i] != 0) {
          carry_cut_list[index2].push_back(index2_child[i]);
          //carry_cut_list[index2_child[i]].push_back(index2);
        }
      }
      //std::cout << index1 << " " << index2 << " does not fit\n";
      return 0;
    }

    // index1 definitely cannot map
    //  and index2 definitely can map
    if (n_unique_index1 > 5) {
      for (uint32_t j = 0; j < n_shared_index2; j++) {
        carry_cut_list[index2].push_back(abde1f1[j]);
        //carry_cut_list[abde1f1[j]].push_back(index2);
      } 
      for (uint32_t i = 0; i < 2; i++) { 
        if (index1_child[i] != 0) {
          carry_cut_list[index1].push_back(index1_child[i]);
          //carry_cut_list[index1_child[i]].push_back(index1);
        }
        if (index2_child[i] != 0) {
          carry_cut_list[index2_child[i]].push_back(index2);
        }
      }
      //std::cout << "Index 1 does not fit\n";
      //std::cout << index1 << " does not fit " << index2 << " does fit\n";
      return 1;
    }

    // index2 definitely cannot map (or is not used)
    //  and index1 definitely can map
    if (index2 == 0 || n_unique_index2 > 5) {
      for (uint32_t j = 0; j < n_shared_index1; j++) {
        carry_cut_list[index1].push_back(abce0f0[j]);
        //carry_cut_list[abce0f0[j]].push_back(index1);
      } 
      for (uint32_t i = 0; i < 2; i++) { 
        // Update carry cut list
        if (index2_child[i] != 0) {
          carry_cut_list[index2].push_back(index2_child[i]);
          //carry_cut_list[index2_child[i]].push_back(index2);
        }
        if (index1_child[i] != 0) {
          carry_cut_list[index1_child[i]].push_back(index1);
        }
      }
      //std::cout << "Index 2 does not fit\n";
      //std::cout << index2 << " does not fit " << index1 << " does fit\n";
      return 1;
    }
  
    // both index have less than 5 unique inputs
    // unique inputs added together - n_matches should be less than 8
    uint32_t n_matches = 0; 
    for (uint32_t i = 0; i < n_shared_index1; i++) {
      for (uint32_t j = 0; j < n_shared_index2; j++) {
        if (abce0f0[i] == abde1f1[j]) n_matches++;
      }
    }

    if (n_unique_index1 + n_unique_index2 - n_matches <= 8) {
      // Update carry cut list
      for (uint32_t j = 0; j < n_shared_index1; j++) {
        carry_cut_list[index1].push_back(abce0f0[j]);
        //carry_cut_list[abce0f0[j]].push_back(index1);
      }
      for (uint32_t j = 0; j < n_shared_index2; j++) {
        carry_cut_list[index2].push_back(abde1f1[j]);
        //carry_cut_list[abde1f1[j]].push_back(index2);
      } 
      for (uint32_t i = 0; i < 2; i++) { 
        if (index1_child[i] != 0) {
          carry_cut_list[index1_child[i]].push_back(index1);
        }
        if (index2_child[i] != 0) {
          carry_cut_list[index2_child[i]].push_back(index2);
        }
      }
      //std::cout << index1 << " " << index2 << " does fit\n";
      return 2;
    }

    for (uint32_t i = 0; i < 2; i++) { 
      if (index1_child[i] != 0) {
        carry_cut_list[index1].push_back(index1_child[i]);
        //carry_cut_list[index1_child[i]].push_back(index1);
      }
      if (index2_child[i] != 0) {
        carry_cut_list[index2].push_back(index2_child[i]);
        //carry_cut_list[index2_child[i]].push_back(index2);
      }
    }
    //std::cout << index1 << " " << index2 << " does not fit\n";

    return 0;
  } 

  std::pair<float, uint32_t> cut_flow( cut_t const& cut )
  {
    uint32_t time{0u};
    float flow{0.0f};

    for ( auto leaf : cut )
    {
      time = std::max( time, delays[leaf] );
      flow += flows[leaf];
    }


    return {flow + cut_area( cut ), time + LUT_DELAY};
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
  void compute_best_cut( uint32_t index )
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

      if ( best_cut == -1 || best_flow > flow + mf_eps || ( best_flow > flow - mf_eps && best_time > time ) )
      {
        best_cut = cut_index;
        best_flow = flow;
        best_time = time;
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
  }

  // Recursively look at the tree
  // Check, is the node one of the cut?
  //    then end there
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

      // In the case of carry node, either map the LUT as
      /*if ( is_a_carry_node(n) ) {
        //if (ps.carry_lut_combined) {
          ntk.add_to_carry_mapping (n, 0); 
          //std::cout << "adding " << n << " to carry mapping.\n";
          for ( auto c: carry_cut_list[ntk.node_to_index(n)]) {
            nodes.push_back(c);
          }
      } else {*/
        if ( map_refs[index] == 0 ) 
          continue;
        for ( auto const& l : cuts.cuts( index ).best() ) {
          nodes.push_back( ntk.index_to_node( l ) );
        }
      //}
      ntk.add_to_mapping( n, nodes.begin(), nodes.end() );

      // This is translated into LUT functionality
      if constexpr ( StoreFunction )
      {
        if ( !is_a_carry_node(n) ) {
          ntk.set_cell_function( n, cuts.truth_table( cuts.cuts( index ).best() ) );
          //std::cout << n << ": ";
          //kitty::print_hex(cuts.truth_table( cuts.cuts( index ).best() ));
          //std::cout << "\n";

          //for (auto cut_node: cuts.cuts(index).best()) {
          //  std::cout << cut_node << " ";
          //} std::cout << "\n";
        }
      }
    }
  }

  void print_state()
  {
    for ( auto i = 0u; i < ntk.size(); ++i )
    {
      std::cout << fmt::format( "*** Obj = {:>3} (node = {:>3})  FlowRefs = {:5.2f}  MapRefs = {:>2}  Flow = {:5.2f}  Delay = {:>3} Carry = {}{}\n", \
        i, ntk.index_to_node( i ), flow_refs[i], map_refs[i], flows[i], delays[i], carry_nodes[i], carry_cut_list[i].empty());
    }
    std::cout << fmt::format( "Level = {}  Area = {}\n", delay, area );



/*    for (auto const& n : top_order) {
      uint32_t i = ntk.node_to_index(n);
      std::cout << fmt::format( "*** Obj = {:>3} (node = {:>3})  FlowRefs = {:5.2f}  MapRefs = {:>2}  Flow = {:5.2f}  Delay = {:>3}\n", i, ntk.index_to_node( i ), flow_refs[i], map_refs[i], flows[i], delays[i] );
    }


    ntk.foreach_po( [this]( auto s ) {
      const auto i = ntk.node_to_index( ntk.get_node( s ) );
      std::cout << fmt::format( "po*** Obj = {:>3} (node = {:>3})  FlowRefs = {:5.2f}  MapRefs = {:>2}  Flow = {:5.2f}  Delay = {:>3}\n", i, ntk.index_to_node( i ), flow_refs[i], map_refs[i], flows[i], delays[i] );
    });*/
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
  std::vector<std::vector<uint32_t>> carry_cut_list;
  std::vector<node<Ntk>> carry_nodes;
  std::vector<node<Ntk>> carry_driver_nodes;
  std::vector<std::vector<int>> carry_cut_index_list;
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
