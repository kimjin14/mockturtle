/*
  post_carry_lut_mapping.hpp

  separate lists for carry nodes and carry lut nodes
  two method to check if node belongs to either carry or carry lut
*/

#pragma once

#include <cstdint>
#include <algorithm>

#include <fmt/format.h>

#include "../utils/stopwatch.hpp"
#include "../views/topo_view.hpp"
#include "../networks/mig.hpp"
#include "cut_enumeration.hpp"
#include "cut_enumeration/mf_cut.hpp"
#include <kitty/print.hpp>

#define MAX_SEARCH 10000000 

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
    cut_enumeration_ps.cut_limit = 8;
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
  bool lut_sharing{false};

  /* Length of paths */
  uint32_t length{5};

  /* Direct path in LUT */
  bool direct{true};

  /* Number of critical path per LUT */
  uint32_t num_crit_path {1};

  /* Delay optimized tech mapping */
  bool delay {false};

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
        carry_paths( ),
        mapped_to_5LUTa( ntk.size(), 0 ),
        mapped_to_5LUTb( ntk.size(), 0 ),
        mapped_to_5LUT( ntk.size(), false),
        mapped_to_5LUT_complemented( ntk.size(), false ),
        num_critical_path_through_node ( ntk.size(), 0 ),
        num_critical_map_refs ( ntk.size(), 0 ),
        direct_input_nodes( ntk.size(), false),
        remapped_nodes( ntk.size() ),
        remapped_nodes_flip( ntk.size(), false ),
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
    print_params();

    set_mapping_refs<false>();

    while ( iteration < ps.rounds )
    {
      compute_mapping<false>(ps.delay, !ps.carry_mapping);
    }
    std::cout << "before ELA\n";
    print_mappable_lut_crit();
    while ( iteration < ps.rounds + ps.rounds_ela)
    {
      compute_mapping<true>(ps.delay, !ps.carry_mapping);
    //set_mapping_refs<true>();
    }
    std::cout << "after ELA\n";
    //print_depth_luts();
    print_mappable_lut_crit();
    //print_mappable_lut_crit(true);
    //print_mappable_lut_crit(true, true);
    
    if (ps.carry_mapping) {
      set_critical_node();
      update_delay();
      uint32_t baseline_delay = delay;
      uint32_t best_delay = delay;
      uint32_t best_area = 0;
      uint32_t best_lut_len = 10;
      uint32_t best_mig_len = 2;
      uint32_t best_max_mig_len = 1000;
      uint32_t best_n_crit = 1;
      uint32_t best_offset = 0;

      //for (uint32_t offset = 0; offset <= LUT_DELAY*10 && offset <= delay; offset+=2*LUT_DELAY) {
      for (uint32_t offset = 0; offset <= 0 && offset <= delay; offset+=LUT_DELAY) {
        for (uint32_t n_crit = 1; n_crit < 7; n_crit++) { // 1 2 
          for (uint32_t mig_len = 2; mig_len <= 10; mig_len*=2) { // 2 4 8 
            for (uint32_t lut_len = 5; lut_len > 1; lut_len/=2) { // 5 2
              for (uint32_t max_mig_len = 100; max_mig_len <= 100; max_mig_len+=20) {// 100 
                std::cout << "Trying config " << offset << " " << lut_len << " " << n_crit << " " << mig_len  << " " << max_mig_len << "\n";   
                map_carry_nodes(false, offset, lut_len, n_crit, mig_len, max_mig_len);
                if ((delay < best_delay) || (delay == best_delay && area < best_area)) {
                  std::cout << "\t-> ";
                  print_state();
                  best_delay = delay;
                  best_area = area;
                  best_lut_len = lut_len;
                  best_mig_len = mig_len;
                  best_max_mig_len = max_mig_len;
                  best_n_crit = n_crit;
                  best_offset = offset;
                } 
                clear_carry_mapping();
              }
            }
          }
        }
      }
      if (best_delay > (int)(baseline_delay*0.95)) std::cout << "\nSelected baseline " << best_delay << " " << baseline_delay << "\n";
      else {
        std::cout << "\nSelected " << best_offset << " " << best_lut_len << " " << best_n_crit << " " << best_mig_len << " " << best_max_mig_len << "\n";   
        map_carry_nodes(true, best_offset, best_lut_len, best_n_crit, best_mig_len, best_max_mig_len);
        std::cout << "There were " << carry_paths.size() << " paths placed.\n";
      }
      //map_carry_nodes(true, 40, 2, 2, 8, 100);
    
      update_delay();
      update_area();
      std::cout << "\t-> ";
      print_state();
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

        init_carry_chain_mapping();
        set_carry_mapping_refs();
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
      init_carry_chain_mapping();
      set_carry_mapping_refs();
      remove_inverter_for_carry_mapping();
      print_state();
      std::cout << "There were " << carry_paths.size() << " paths placed on carry chain\n"; 
  }

  ///////////////////////////////////////////////////////////////
  // Selecting Path to be Placed on Carry Chain and Mapping
  ///////////////////////////////////////////////////////////////

  void count_mappable_lut_crit( uint32_t index, uint32_t& total_num_cuts, uint32_t& total_num_mappable_cuts, bool remap = false , bool direct = false) {


    auto n = ntk.index_to_node(index); 
    if ( ntk.is_constant( n ) || ntk.is_pi( n ) || ntk.visited( n ) == 1 ) {
      return;
    }

    //std::cout << "index " << index << "(" << delays[index] << ")\n";
    total_num_cuts++;
    ntk.set_visited(n, 1);

    uint32_t count_crit = 0;
    uint32_t count_mappable_crit = 0;

    for ( auto leaf: cuts.cuts(index)[0] ) {
      if (delays[leaf] == (delays[index] - LUT_DELAY)) {
        
        //if ( count_path_to_node (index, index, leaf, 0, false) == 1  || (remap && try_other_graph(index, leaf)) || (direct && check_direct_input(leaf,index)) ) {
            //&& length_to_node (index, index, leaf) < 1000) { //direct
            //check_direct_input (leaf, index) || 
            //check_shorter_path(leaf, index, index )) {`
            count_mappable_crit++; 
        //}
        count_mappable_lut_crit(leaf, total_num_cuts, total_num_mappable_cuts, remap, direct);
        count_crit++;
      }
    }
    if (count_crit < 2 && count_mappable_crit > 0 ) total_num_mappable_cuts++;
  }

  void print_mappable_lut_crit ( bool remap = false, bool direct = false ) {

    uint32_t total_num_cuts = 0;
    uint32_t total_num_mappable_cuts = 0;
    ntk.clear_visited();
    ntk.foreach_po( [&]( auto s ) {

      uint32_t n = ntk.get_node(s);
      uint32_t index = ntk.node_to_index(n);
  
      if ( ntk.is_constant( n ) || ntk.is_pi( n ) ) {
        return;
      }

      if (delay == delays[index]) {
        count_mappable_lut_crit (index, total_num_cuts, total_num_mappable_cuts, remap, direct);
      }
    });
    std::cout << "Mappable nodes " << total_num_mappable_cuts << " out of " << total_num_cuts << ".\n";
  }

  void print_depth_luts () {
    uint32_t n_out_accum = 0;
    uint32_t n_out = 0;
    uint32_t n_worst = 0;
    for ( uint32_t d = 0; d < delay+LUT_DELAY; d+=LUT_DELAY)
    {
      n_out = 0;
      ntk.foreach_po( [&]( auto s ) {

        uint32_t n = ntk.get_node(s);
        uint32_t index = ntk.node_to_index(n);
  
       // if ( ntk.is_constant( n ) || ntk.is_pi( n ) ) {
       //   return;
       // }
        
        if ( d == delays[index] ) {
          if (d == delay) n_worst++;
          n_out++;
        }
      });
      n_out_accum += n_out;

      //std::cout << "Level =\t" << d << ".\t";
      //std::cout << "COs =\t" << n_out << ".\t";
      //std::cout << "\t" << std::fixed << std::setprecision(0) << (float)n_out_accum/ntk.num_pos()*100 << "\%\n";
      std::cout << std::fixed << std::setprecision(0) << (float)n_out_accum/ntk.num_pos()*100 << "\n";
    } 

    std::cout << "Worst LUT level " << delay << ", COs = " << n_worst << "/" << ntk.num_pos() << "\n";
 
  }

  uint32_t lut_mappable ( uint32_t index ) {

    uint32_t n_crit = 0;
    uint32_t n_mappable_crit = 0;

    for (auto leaf: cuts.cuts(index)[0]) {
      if (delays[leaf] == (delays[index] - LUT_DELAY)) {
        if (count_path_to_node(index, index, leaf, 0) == 1) n_mappable_crit++;
        n_crit++; 
      }
    }

    //if (n_crit > 2) return false;
    if ((n_crit - n_mappable_crit) == 0) return true;
    return false;
  }

  uint32_t length_to_node (  uint32_t index, uint32_t source_index,  uint32_t dest_index) {

    auto const& n = ntk.index_to_node(source_index);

    if (source_index == dest_index) {
      return 1;
    }

    /*if (is_a_carry_node(ntk.index_to_node(index))) {
      for ( auto leaf : carry_cut_list[index] )
        if (source_index == leaf) return 0;
    } else {*/
      for ( auto leaf : cuts.cuts(index)[0] )
        if (source_index == leaf) return 0;
    //}
 
    uint32_t total = 0;

    for (uint32_t i = 0; i < ntk.fanin_size(n); i++) {
      auto const leaf = ntk.get_children(n,i);
      auto const leaf_index = ntk.node_to_index(leaf);
      //if (!is_a_carry_node(leaf) && !ntk.is_constant(leaf)) {
      if (!ntk.is_constant(leaf)) {
        //std::cout << source_index << " " << leaf_index << " " << dest_index << "\n";
        total = length_to_node(index, leaf_index, dest_index);
        if (total > 0) break;
      }
    }
    if ( total > 0 ) return total+1;
    return 0;


  }
  // returns true if no nodes in the path is on the carry 
  bool no_carry_on_path(  uint32_t index, uint32_t source_index,  uint32_t dest_index, uint32_t use_remap = false) {

    auto  n = ntk.index_to_node(source_index);

    if(use_remap && source_index < tmp_size && remapped_nodes[index][source_index] > -1) {
      std::cout << "from " << n << "(" << remapped_nodes[index][source_index] << ")";
      n = ntk.index_to_node(ntk.get_node(tmp_signals[remapped_nodes[index][source_index]]));
      source_index = ntk.node_to_index(n);
      std::cout << " to " << n;
      std::cout << "\n";
    }

    if (source_index == dest_index) {
      return true;
    } 

    if (is_a_carry_node(ntk.index_to_node(dest_index))) return false;
    
    if (is_a_carry_node(ntk.index_to_node(index))) {
      for ( auto leaf : carry_cut_list[index] )
        if (source_index == leaf) return false;
    } else {
      for ( auto leaf : cuts.cuts(index)[0] )
        if (source_index == leaf) return false;
    }
 
    for (uint32_t i = 0; i < ntk.fanin_size(n); i++) {
      auto const leaf = ntk.get_children(n,i);
      auto const leaf_index = ntk.node_to_index(leaf);
      if (!ntk.is_constant(leaf))
        if (no_carry_on_path(index, leaf_index, dest_index, use_remap)){
          //std::cout << "\t" << source_index << " is " << is_a_carry_node(source_index) << "\n";
          return !is_a_carry_node(source_index);
        }
    }
    return false;
  }


  uint32_t count_path_to_node (  uint32_t index, uint32_t source_index,  uint32_t dest_index, uint32_t cut_index, bool use_remap = false) {

    auto n = ntk.index_to_node(source_index);

    if(use_remap && source_index < tmp_size && remapped_nodes[index][source_index] > -1) {
      //std::cout << "from " << n << "(" << remapped_nodes[index][source_index] << ")";
      n = ntk.index_to_node(ntk.get_node(tmp_signals[remapped_nodes[index][source_index]]));
      source_index = ntk.node_to_index(n);
      //std::cout << " to " << n;
      //std::cout << "\n";
    }

    if (source_index == dest_index) {
      return 1;
    }

    /*if (is_a_carry_node(ntk.index_to_node(index))) {
      for ( auto leaf : carry_cut_list[index] )
        if (source_index == leaf) return 0;
    } else {*/
      for ( auto leaf : cuts.cuts(index)[cut_index] )
        if (source_index == leaf) return 0;
    //}
 
    uint32_t total = 0;

    for (uint32_t i = 0; i < ntk.fanin_size(n); i++) {
      auto const leaf = ntk.get_children(n,i);
      auto const leaf_index = ntk.node_to_index(leaf);
      if (!ntk.is_constant(leaf))
        total += count_path_to_node(index, leaf_index, dest_index, cut_index, use_remap);
    }
    return total;
  }


  // Depth first search 
  // Find path from output of LUT to target input LUT
  bool add_path_within_LUT_to_carry (std::vector<node<Ntk>>& path_for_carry_chain, 
      uint32_t index, uint32_t source_index, uint32_t dest_index) {

    auto n = ntk.index_to_node(source_index);
    uint32_t orig_index = source_index;

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

    if( source_index < tmp_size && remapped_nodes[index][source_index] > -1 ) {
      //    remapped_nodes[index][index] == -2 && source_index != index ) {
      //std::cout << "from " << source_index << "(" << remapped_nodes[index][source_index] << ")";
      n = ntk.get_node(tmp_signals[remapped_nodes[index][source_index]]);
      source_index = ntk.node_to_index(n);
      //std::cout << " to " << source_index;
      //std::cout << "\n";
    } 
 
    bool found = false; 
    for (uint32_t i = 0; i < ntk.fanin_size(n); i++) {
      auto const leaf = ntk.get_children(n,i);
      auto const leaf_index = ntk.node_to_index(leaf);
      //if (!is_a_carry_node(leaf) && !ntk.is_constant(leaf))
      if (!ntk.is_constant(leaf))
        found = add_path_within_LUT_to_carry (path_for_carry_chain, index, leaf_index, dest_index);
      //else std::cout << leaf << " carry or const?\n"; 

      if (found) {
        if (orig_index == index) path_for_carry_chain.push_back(orig_index);
        else path_for_carry_chain.push_back(source_index);
        if (ps.verbose && ps.verbosity > 3) std::cout << source_index <<  " -> " << dest_index << " ";
        return true;
      }
      //else std::cout << leaf << " not found?\n"; 
    }
    return found;
  }
  bool add_cut_to_carry(std::vector<node<Ntk>>& path_for_carry_chain, 
      uint32_t index, uint32_t source_index, uint32_t dest_index) {

    auto n = ntk.index_to_node(source_index);
    uint32_t orig_index = source_index;

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
    if( source_index < tmp_size && remapped_nodes[index][source_index] > -1) {//||
      //    remapped_nodes[index][index] == -2 && source_index != index ) {
      //std::cout << "from " << curr_index << "(" << remapped_nodes[index][curr_index] << ")";
      n = ntk.get_node(tmp_signals[remapped_nodes[index][source_index]]);
      source_index = ntk.node_to_index(n);
      //std::cout << " to " << curr_index;
      //std::cout << "\n";
    }

    bool found = false; 
    for (uint32_t i = 0; i < ntk.fanin_size(n); i++) {
      auto const leaf = ntk.get_children(n,i);
      auto const leaf_index = ntk.node_to_index(leaf);
      if (!ntk.is_constant(leaf))
        found = add_cut_to_carry(path_for_carry_chain, index, leaf_index, dest_index);
      if (found) {
        if(orig_index == index) {
          add_dependent_cut_to_carry(index, orig_index, source_index, leaf_index);
          carry_cut_list[orig_index].push_back(leaf_index);
        } else  {
          add_dependent_cut_to_carry(index, source_index, source_index, leaf_index);
          carry_cut_list[source_index].push_back(leaf_index);
        }
        return true;
      }
    }
    return found;
  }

  // Adds leafs that are part of the node
  // index: used to find relevant cut
  // source_index: adding the subset of cut to here
  // curr_index: index used to get to leaves
  void add_dependent_cut_to_carry(uint32_t index, uint32_t source_index, uint32_t curr_index, uint32_t carry_index) {
  

    for (auto leaf: cuts.cuts(index)[0]) {
      if (leaf == curr_index) {
        bool in_list = false;
        for (auto inlist: carry_cut_list[source_index]) {
          if (curr_index == inlist) in_list = true;
        }
        if(!in_list && curr_index) carry_cut_list[source_index].push_back(curr_index);
        return;
      }
    }

    for (uint32_t i = 0; i < ntk.fanin_size(curr_index); i++) {
      auto const child = ntk.get_children(ntk.index_to_node(curr_index),i);
      if (!ntk.is_constant(child) && ntk.node_to_index(child) != carry_index)
        add_dependent_cut_to_carry(index, source_index, ntk.node_to_index(child), carry_index); 
    }

    return;
  }

  uint32_t select_next_node (uint32_t index, auto cut_list, bool useeasy) {

    uint32_t max_leaf = index;
    uint32_t max_delay = 0;

    //uint32_t mult_max_delay = 0;
    for ( auto leaf : cut_list ) {
      if (!is_a_carry_node(leaf)) { 
        if (ps.verbose && ps.verbosity > 2) std::cout << "\t" << leaf << ":" << delays[leaf] << "," << flow_refs[leaf] << "," << map_refs[leaf] <<"\n"; 
        if (max_delay < delays[leaf]) {
          max_leaf = leaf;
          max_delay = delays[leaf];
        }
      } else 
        if (ps.verbose && ps.verbosity > 2) std::cout << "\tx" << leaf << ":" << delays[leaf] << "," << flow_refs[leaf] << "," << map_refs[leaf] <<"\n"; 
    }

    // Check that the input feeds the output directly (is one of the child node)
    // And it doesn't feed anything else
    bool easy = false;
    ntk.foreach_fanin( ntk.node_to_index(index), [&]( auto const& f ) {
      if (ntk.get_node(f) == max_leaf) {
        uint32_t total = count_path_to_node (index, index, max_leaf, 0);
        std::cout << "\t\tcount is " << total << "\n";
        if (total == 1)
          easy = true; 
      }
    });

    if (useeasy && !easy) {
      //max_leaf = index;
      std::cout << "\t\tNOT EASY\n";
    }
    return max_leaf;
  }

  // Keep going deeper to its fanin until depth is found
  bool find_deepest_LUT (std::vector<node<Ntk>>& path_for_carry_chain, uint32_t index, uint32_t length) { 
    
    if (ps.verbose && ps.verbosity > 2) std::cout << "Index " << index << "(" << delays[index] << "):";
    //if (delays[index] == 20 &&  length > 6) {
    if (delays[index] == 20 && length >= 5) {
      auto leaf = cuts.cuts( index )[0].begin()[0];
      if(!add_path_within_LUT_to_carry(path_for_carry_chain, index, index, leaf)) {
        if (ps.verbose && ps.verbosity > 2) std::cout << "\n\t\tnot added\n";
        return false;
      } 
      return true;
    }

    bool deepest = false;

    uint32_t max_leaf = index;

    if (is_a_carry_node(ntk.index_to_node(index))) {
      if (ps.verbose && ps.verbosity > 2) std::cout << " carry\n"; 
      max_leaf = select_next_node(index, carry_cut_list[index], false);
      if (max_leaf != index) {
        deepest = find_deepest_LUT(path_for_carry_chain, max_leaf, 1);

        if (deepest) {
          path_for_carry_chain.push_back(max_leaf);
        }
      }
    } else {
      if (ps.verbose && ps.verbosity > 2) std::cout << " lut\n"; 
      max_leaf = select_next_node(index, cuts.cuts( index )[0], false );
      if (max_leaf != index) {
        deepest = find_deepest_LUT(path_for_carry_chain, max_leaf, length+1);

        if (deepest)  {
          std::cout << "count is " << count_path_to_node (index, index, max_leaf, 0) << "\n";
          ntk.clear_visited();
          if(add_path_within_LUT_to_carry(path_for_carry_chain, index, index, max_leaf)) {
            if (ps.verbose && ps.verbosity > 2) std::cout << "\t\tAdded path\n";
          } else { //cannot put this in carry
            return false;
          } 
        }
      } 
    }

    return deepest;
  }

  // Keep going deeper to its fanin until depth is found
  bool find_easiest_LUT (uint32_t& counter, std::vector<node<Ntk>>& path_for_carry_chain, uint32_t index, \
      uint32_t length, uint32_t mig_length, uint32_t& pi_count, uint32_t& chain_len_count, uint32_t& crit_count, uint32_t& mig_len_count,\
       uint32_t& mult_path_count, uint32_t lut_chain_len, uint32_t num_crit_lut_in, uint32_t mig_len, uint32_t max_mig_len) { 

    //if (counter > MAX_SEARCH || ntk.visited(ntk.index_to_node(index)) > 0) {
    //  if (ps.verbose && ps.verbosity > 3) std::cout << "\t\t\tCOUNTER/VISITED " << index << " " << length << "\n";
    //  return false; 
    //}

    if (ps.verbose && ps.verbosity > 3) 
      std::cout << "\t\tIndex " << index << "(" << delays[index] << "): " << length << "\n";

    counter++;   

    if ( mig_length >= max_mig_len ) {
      //std::cout << "\t\tMIG LENGTH IS TOO LONG " << index << "(" << delays[index] << "): " << mig_length << "\n";
      path_for_carry_chain.push_back(index);
      return true;
    }
  
    // If first LUT, add any input 
    if ( counter > MAX_SEARCH || ntk.visited(ntk.index_to_node(index)) > 0 || delays[index] == 20 ) {
      pi_count++;
      if (ps.verbose && ps.verbosity > 3) std::cout << "\t\t\tPI " << index << " " << length << "\n";
      if (length >= lut_chain_len) {
        path_for_carry_chain.push_back(index);
        return true;
      }
      chain_len_count++;
      return false; 
    }
    
    std::vector<uint32_t>  cut_list(cuts.cuts(index)[0].begin(), cuts.cuts(index)[0].end());
    if (is_a_carry_node(index)) {
      cut_list = carry_cut_list[index];
    } 

    // Check 1: checks that it has only 1 critical input
    uint32_t crit_child_count = 0;     
    uint32_t crit_child_node = 0;
    for ( auto leaf: cut_list ) { 
      if ((delays[index] - LUT_DELAY) == delays[leaf] || (delays[index] - CARRY_DELAY) == delays[leaf] || (delays[index] - LUT_ADDER_DELAY) == delays[leaf]) {
        crit_child_count++;
        crit_child_node = leaf;
      }
    }
    if (crit_child_count > num_crit_lut_in) {
      ntk.set_visited(ntk.index_to_node(index),1);
      crit_count++;
      if (ps.verbose && ps.verbosity > 3) std::cout << "\t\t\tMULT CRIT " << crit_child_node << " " << crit_child_count << " " << length <<"\n";
      if (length >= lut_chain_len) {
        path_for_carry_chain.push_back(index);
        return true;
      } 
      chain_len_count++;
      return false;
    }
    for ( auto leaf: cut_list ) { 
      if (((delays[index] - LUT_DELAY) == delays[leaf] || (delays[index] - CARRY_DELAY) == delays[leaf] || (delays[index] - LUT_ADDER_DELAY) == delays[leaf]) && path_for_carry_chain.empty()) {
        crit_child_node = leaf;

        // Check 1: only one path to the leaf node  
        //uint32_t num_path = count_path_to_node (index, index, crit_child_node, 0);
        //uint32_t num_path = 1;
        // Check 2: not on carry already (all nodes in that particular path)
        //bool on_carry = !no_carry_on_path(index, index, crit_child_node, true);
        bool on_carry = is_a_carry_node(index) || is_a_carry_node(crit_child_node);
        // Check 3: mig length within LUT
        //uint32_t len_path = length_to_node (index, index, crit_child_node);
        //uint32_t len_path = 1;
        //bool path_too_long = len_path > mig_len; 

        //bool direct_input = !on_carry && check_direct_input (crit_child_node, index);
        //bool shorter_input = !on_carry && check_shorter_path (crit_child_node, index, index);
        //if (num_path > 1 && !direct_input && !on_carry && !path_too_long && try_other_graph(index, crit_child_node)) {
        //  num_path = 1;
        //}
        //bool mult_path_through_lut = num_path != 1;
        // Required for legal mapping
        if ( on_carry ) {
          //if ( ps.verbose && ps.verbosity > 3 && on_carry ) std::cout << "\t\t\tCANNOT BE PLACED " << crit_child_node << " is on carry\n";
          //else if ( ps.verbose && ps.verbosity > 3 ) std::cout << "\t\t\tCANNOT BE PLACED " << crit_child_node << " shorter?" <<  shorter_input << " n_path?" << num_path << " path_len?" << len_path << " mig_len?" << mig_len << "\n";

          if (length >= lut_chain_len) {
            path_for_carry_chain.push_back(index);
            return true;
          } 
          chain_len_count++;
          if (find_easiest_LUT(counter, path_for_carry_chain, crit_child_node, 0, 0,\
              pi_count, chain_len_count, crit_count, mig_len_count, mult_path_count,
              lut_chain_len, num_crit_lut_in, mig_len, max_mig_len))
            return false;
          else continue;
        }
        if (find_easiest_LUT(counter, path_for_carry_chain, crit_child_node, length+1, mig_length+mig_len,\
            pi_count, chain_len_count, crit_count, mig_len_count, mult_path_count,
            lut_chain_len, num_crit_lut_in, mig_len, max_mig_len)) {
          // add node
          path_for_carry_chain.push_back(index);
          // add cut
          carry_cut_list[index].clear();
          for (auto icut: cuts.cuts(index)[0]) {
            if (icut != crit_child_node) {
              cut_list.push_back(icut);
              carry_cut_list[index].push_back(icut);
            }
          } 
          carry_cut_list[index].push_back(crit_child_node);
          direct_input_nodes[index] = true;
          return true;
        } else {
          if (!is_a_carry_node(index))remapped_nodes[index][index] = -1;
        } 
      }
    }
    ntk.set_visited(ntk.index_to_node(index),1);

    return false;
  } 
  void print_lut (uint32_t index, uint32_t curr_index, bool use_remap = false) {

    if (ntk.is_constant(curr_index)) return;
    for (auto leaf: cuts.cuts(index)[0]) {
      if (curr_index == leaf) return;
    }

    if(use_remap && curr_index < tmp_size && remapped_nodes[index][curr_index] > -1 ) {
      //std::cout << "from " << curr_index << "(" << remapped_nodes[index][curr_index] << ")";
      auto n = ntk.get_node(tmp_signals[remapped_nodes[index][curr_index]]);
      curr_index = ntk.node_to_index(n);
      //std::cout << " to " << curr_index;
      //std::cout << "\n";
    }
    std::cout << curr_index << "->";
    
    for (uint32_t i = 0; i < ntk.fanin_size(curr_index); i++) {
      auto leaf = ntk.get_children(curr_index, i); 
      bool flipped = false;
      if(use_remap && leaf < tmp_size && remapped_nodes[index][leaf] > -1) {
        leaf = ntk.node_to_index(ntk.get_node(tmp_signals[remapped_nodes[index][leaf]]));
        flipped = true;
      }
      if (ntk.is_complemented_children(curr_index, i)&&!flipped) std::cout << "*"<< leaf << " "; 
      else if (!ntk.is_complemented_children(curr_index, i)&&flipped) std::cout << "*"<< leaf << " "; 
      else std::cout << leaf << " "; 
    }    
    std::cout << "\n";
 
    for (uint32_t i = 0; i < ntk.fanin_size(curr_index); i++) {
      auto leaf = ntk.get_children(curr_index, i); 
      print_lut (index, leaf, use_remap);
    }
    return;
  }
/*  void print_lut (uint32_t index, uint32_t curr_index) {


    for (auto leaf: cuts.cuts(index)[0]) {
      if (curr_index == leaf) return;
    }

    for (uint32_t i = 0; i < ntk.fanin_size(curr_index); i++) {
      auto leaf = ntk.get_children(curr_index, i); 
      std::cout << curr_index << "->" << leaf << "\n";    
      if (ntk.is_constant(leaf)) continue;
      print_lut (index, leaf);
    }
    return;
  }*/

  void map_carry_nodes (bool set, uint32_t offset, uint32_t lut_chain_len, uint32_t num_crit_lut_in, uint32_t mig_len, uint32_t max_mig_len) {

    uint32_t pi_count = 0;
    uint32_t chain_len_count = 0;
    uint32_t crit_count = 0;
    uint32_t mig_len_count = 0;
    uint32_t mult_path_count = 0;
    uint32_t counter = 0;

    ntk.clear_visited();

    // Search for POs with delays that are within 2 LUT delay of critical delay
    //for (uint32_t delay_offset = 0; delay_offset <= 2*LUT_DELAY && counter < MAX_SEARCH;){
    for (uint32_t delay_offset = 0; delay_offset <=  offset&& counter < MAX_SEARCH;){
      if (ps.verbose && ps.verbosity > 3) std::cout << "\tMapping " << delay << " to " << delay-delay_offset << "\n";
      bool mapped_something = false;
      std::vector<node<Ntk>> path_for_carry_chain;
      ntk.foreach_po( [&]( auto const& s ) {
        const auto node = ntk.get_node(s);
        const auto index = ntk.node_to_index(node);
         
        if (ntk.is_pi(node) || ntk.is_constant(node) || counter > MAX_SEARCH) return;

        if (delays[index] >= delay-delay_offset && delays[index] <= delay ) {
          if (ps.verbose && ps.verbosity > 3) 
            std::cout << "\t\tPO " << index << "(" << delays[index] << "): " << counter << "\n";

          //for(uint32_t i = 0; i < ntk.fanin_size(index); i++) {
          //  auto const child = ntk.get_children(node,i);
            find_easiest_LUT(counter, path_for_carry_chain, index/*ntk.node_to_index(child)*/, 0, 0, \
              pi_count, chain_len_count, crit_count, mig_len_count, mult_path_count, \
              lut_chain_len, num_crit_lut_in, mig_len, max_mig_len);

            if (path_for_carry_chain.empty() || counter > MAX_SEARCH) return;// continue;
            mapped_something = true;
            carry_paths.push_back(path_for_carry_chain);
            for (uint32_t j = 1; j < path_for_carry_chain.size(); j++) {
              carry_nodes[path_for_carry_chain[j]] += 1;
              carry_driver_nodes[path_for_carry_chain[j]] = path_for_carry_chain[j-1];
            }
            path_for_carry_chain.clear();
            update_delay();
            update_area(); 
            ntk.clear_visited();
          if (ps.verbose && ps.verbosity > 3){ 
            print_critical_path();
            std::cout << "\tMapped one chain, updated delay," << delay << "," << area << "\n\n";
          }
        }

      } );
      if (!mapped_something)
        delay_offset+=LUT_DELAY;
    }
    update_delay();
    update_area(); 
    
    if (set) {
      ntk.resize_mapping(tmp_signals.size());
      remove_inverter();
      for (auto carry_path: carry_paths) {
        for (int i = carry_path.size()-1; i >= 1; i--) {
          //if (carry_cut_list[carry_path[i]].empty()) set_direct_input(carry_path[i-1], carry_path[i]);
          if (direct_input_nodes[carry_path[i]]) set_fake_direct_input(carry_path[i-1], carry_path[i]);
          //if (direct_input_nodes[carry_path[i]]) set_direct_input(carry_path[i-1], carry_path[i]);
          //else set_tt_easy_node(carry_path[i-1], carry_path[i]);
        }
      }
      check_inverter();
    }

  }

  // Keep going deeper to its fanin until depth is found
  bool find_all_LUT (uint32_t& counter, std::vector<node<Ntk>>& path_for_carry_chain, uint32_t index, \
      uint32_t length, uint32_t mig_length, uint32_t& pi_count, uint32_t& chain_len_count, uint32_t& crit_count, uint32_t& mig_len_count,\
       uint32_t& mult_path_count, uint32_t lut_chain_len, uint32_t num_crit_lut_in, uint32_t mig_len) { 

    if (counter > MAX_SEARCH || ntk.visited((ntk.index_to_node(index)) > 0)) return false; 

    if (ps.verbose && ps.verbosity > 3) 
      std::cout << "\t\tIndex " << index << "(" << delays[index] << "): " << length << "\n";

    counter++;   


    // If first LUT, add any input 
    if ( ntk.is_pi(index) ) {
      if (ps.verbose && ps.verbosity > 3) std::cout << "\t\t\tPI " << index << " " << length << "\n";
      if (length >= lut_chain_len) {
        path_for_carry_chain.push_back(index);
        return true;
      }
      return false; 
    }
   

    uint32_t worst_delay = 0;
    uint32_t worst_delay_index = 0;
    uint32_t crit_child_count = 0;
    uint32_t crit_child_index = 0;
    uint32_t crit_cost = 0;

    std::cout << "Considering";
    if (is_a_carry_node(index)) {
      std::cout << "x";
      for (uint32_t i = 0; i < carry_cut_list[index].size()-1; i++) {
        auto leaf = carry_cut_list[index][i];
        std::cout << " " << leaf;
        if (delays[leaf] > worst_delay) {
          worst_delay = delays[leaf];
          worst_delay_index = leaf;
        }
      }
    } else {
      for ( auto leaf: cuts.cuts(index).best() ) { 
        std::cout << " " << leaf;
        if (delays[leaf] > worst_delay) {
          worst_delay = delays[leaf];
          worst_delay_index = leaf;
        }
        bool mult_path = count_path_to_node (index, index, leaf, 0) != 1;
        bool on_carry = !no_carry_on_path(index, index, leaf);
        if ( mult_path || on_carry || ntk.is_constant(leaf) || ntk.is_pi(leaf) ) {
          std::cout << "x";
          continue;
        }

        uint32_t cost = delays[leaf];
        //uint32_t cost = map_refs[leaf];
        if ( crit_child_count == 0 || delays[leaf] > delays[crit_child_index] || (delays[leaf] == delays[crit_child_index] && map_refs[leaf] < map_refs[crit_child_index])) {
          crit_child_index = leaf;
          std::cout << "*";
          crit_cost = cost;
        }
        crit_child_count++;
      }  
    }
    std::cout << "\n";

    
    // End if no possible leaf
    if ( crit_child_count == 0 || is_a_carry_node(index) ) {
      //std::cout << " -- w" << worst_delay_index << "\n";
      if (ps.verbose && ps.verbosity > 3) std::cout << "\t\t\tNO POSSIBLE LEAF " << length << "\n";
      if (length >= lut_chain_len) {
        path_for_carry_chain.push_back(index);
        return true;
      }
      if (ntk.is_pi(worst_delay_index) || ntk.is_constant(worst_delay_index)) return false;
      find_all_LUT(counter, path_for_carry_chain, worst_delay_index, 0, 0,\
             pi_count, chain_len_count, crit_count, mig_len_count, mult_path_count,
             lut_chain_len, num_crit_lut_in, mig_len);
      return false; 
    }

    // Place or continue (if carry) 
    //std::cout << " -- c" << crit_child_index << "\n";
    assert(crit_child_index != 0);
    if (find_all_LUT(counter, path_for_carry_chain, crit_child_index, length+1, mig_length,\
        pi_count, chain_len_count, crit_count, mig_len_count, mult_path_count,
        lut_chain_len, num_crit_lut_in, mig_len)) {
      if (!add_path_within_LUT_to_carry(path_for_carry_chain, index, index, crit_child_index) ){
        std::cout << "DID NOT PLACE ANYTHING WHY?\n";
        assert(0);
      }
      add_cut_to_carry(path_for_carry_chain, index, index, crit_child_index);  
      return true;
    } 
    //ntk.set_visited(ntk.index_to_node(index),1);

    return false;
  } 

  void map_all_carry_nodes (bool set, uint32_t lut_chain_len, uint32_t num_crit_lut_in, uint32_t mig_len) {

    uint32_t old_delay = delay;
    uint32_t pi_count = 0;
    uint32_t chain_len_count = 0;
    uint32_t crit_count = 0;
    uint32_t mig_len_count = 0;
    uint32_t mult_path_count = 0;
    uint32_t counter = 0;

    ntk.clear_visited();

    std::vector<node<Ntk>> path_for_carry_chain;

    for ( uint32_t delay_offset = 0; delay_offset < delay; delay_offset += LUT_DELAY ) {
      if (ps.verbose && ps.verbosity > 3) std::cout << "\tMapping " << delay << " to " << delay-delay_offset << "\n";
      counter = 0;
      ntk.foreach_po( [&]( auto const& s ) {
        const auto node = ntk.get_node(s);
        const auto index = ntk.node_to_index(node);
         
        if (ntk.is_pi(node) || ntk.is_constant(node) || counter > MAX_SEARCH) return;

        if (delays[index] >= delay-delay_offset && delays[index] <= delay ) {

          if (ps.verbose && ps.verbosity > 3) 
            std::cout << "\t\tPO " << index << "(" << delays[index] << "): " << counter << "\n";

          do {
            path_for_carry_chain.clear();
            ntk.clear_visited();
            find_all_LUT(counter, path_for_carry_chain, index, 0, 0, \
                pi_count, chain_len_count, crit_count, mig_len_count, mult_path_count, \
                lut_chain_len, num_crit_lut_in, mig_len);

            if (path_for_carry_chain.empty() || counter > MAX_SEARCH) return;
            carry_paths.push_back(path_for_carry_chain);
            for (uint32_t j = 1; j < path_for_carry_chain.size(); j++) {
              carry_nodes[path_for_carry_chain[j]] += 1;
              carry_driver_nodes[path_for_carry_chain[j]] = path_for_carry_chain[j-1];
            }
            update_delay();
            update_area(); 
      
            if (ps.verbose && ps.verbosity > 3) 
              std::cout << "\tMapped one chain, updated delay," << delay << "," << area << "\n\n";
          } while (!path_for_carry_chain.empty());
        }
      } );
    }

    print_carry_paths();
    update_delay();
    update_area();
    std::cout << "Updated delay is " << delay << "\n";
    std::cout << "Total nodes explored is " << counter << " and failed for " \
      << pi_count << " " \
      << chain_len_count << " " \
      << crit_count << " " \
      << mig_len_count << " " \
      << mult_path_count << "\n";

    
    /*if (old_delay + LUT_DELAY <= delay) {
      std::cout << "WORSE MAPPING CLEAR\n";
      clear_carry_mapping();
    } */
    
    if (set) {
      remove_inverter();
      for (auto carry_path: carry_paths) {
        for (int i = carry_path.size()-1; i >= 1; i--) {
          set_tt_easy_node(carry_path[i-1], carry_path[i]);
        }
      }
      check_inverter();
    }

  }


  void clear_carry_mapping() {
    for (auto carry_path: carry_paths) {
      for (auto cnode: carry_path) {
        carry_nodes[cnode] = 0;
        carry_cut_list[cnode].clear();
        carry_driver_nodes[cnode] = 0;
      }
      carry_path.clear();
    }
    carry_paths.clear();
    update_delay();

 
    for (uint32_t i = 0; i < tmp_size; i++) {
      for (uint32_t j = 0; j < tmp_size; j++) {
        remapped_nodes[i][j] = -1;
      }
      direct_input_nodes[i] = false;
      remapped_nodes_flip[i] = false;
    }
    tmp_signals.clear();
  }
 
 
  void set_tt_easy_node (uint32_t cindex, uint32_t index) {
    //auto cut_leaf_list = cuts.cuts( LUTindex )[0];
    uint32_t cut1_index[5] = {0, 1, 2, 3, 4};
    uint32_t cut2_index[5] = {0, 1, 2, 3, 4};
    if (ps.verbose && ps.verbosity > 2) {
      std::cout << index << ":";
      for (auto icut: carry_cut_list[index]) {
        if (icut != cindex) {
        }
        std::cout << icut << " ";
      }
      std::cout << "\n";
    }
    std::vector<uint32_t> cut_list;
    for (auto icut: carry_cut_list[index]) {
      if (icut != cindex) {
        cut_list.push_back(icut);
      }
    }
    // Get relevant children node for mapping
    node<Ntk> i_child[2] = {0};

    auto tmp_index = index; 
    bool i_child_flipped[2] = {0};
    if (index < tmp_size && remapped_nodes[index][index] > -1) { 
      tmp_index = ntk.node_to_index(ntk.get_node(tmp_signals[remapped_nodes[index][index]]));
    }

    get_children_node (i_child, tmp_index, cindex); 

    //std::cout << "\t" << i_child[0] << " ";
    auto function1 = flip_carry_user_tt(index, i_child[0], cut_list, true, remapped_nodes_flip[index]);
    //kitty::print_hex (function1);
    //std::cout << "\n";
    //std::cout << "\t" << i_child[1] << " ";
    auto function2 = flip_carry_user_tt(index, i_child[1], cut_list, true, remapped_nodes_flip[index]);
    //kitty::print_hex (function2);
    //std::cout << "\n";


    bool child_complement[3] = {0};
    determine_child_complement (child_complement, ntk.index_to_node(tmp_index),\
        cindex, i_child[0], i_child[1]);

    uint32_t isize = carry_cut_list[index].size()-1;
    kitty::dynamic_truth_table function = combine_carry_LUT_function(isize, isize, isize, function1, function2, cut1_index, cut2_index, child_complement);

    if (ps.verbose && ps.verbosity > 2) {
      std::cout << "\tinv:" << child_complement[0] << " " << child_complement[1] << " " << child_complement[2] << "\n";
      std::cout << "\tfinal:";
      kitty::print_hex(function);
      std::cout << "\n";
    }

    ntk.set_cell_function(ntk.index_to_node(tmp_index), function);

  }
  // set correct mask if possible
  //currently not correct
   bool set_fake_direct_input (uint32_t cindex, uint32_t index) {
    std::vector<uint32_t> cut_list;
    for (auto icut: cuts.cuts(index)[0]) {
      if (icut != cindex) {
        cut_list.push_back(icut);
      }
    } 

    auto fcn_xsub1 = set_tt_const(index, ntk.get_children(index, 0), cut_list, cindex, 0);
    auto fcn_xsub2 = set_tt_const(index, ntk.get_children(index, 1), cut_list, cindex, 0);
    auto fcn_sub1 = ntk.is_complemented_children(index,0) ? ~fcn_xsub1:fcn_xsub1; 
    auto fcn_sub2 = ntk.is_complemented_children(index,1) ? ~fcn_xsub2:fcn_xsub2;
    auto fcn_sub_OR1 = kitty::binary_or ( fcn_sub1, fcn_sub2 ); 
    auto fcn_sub_AND1 = kitty::binary_and ( fcn_sub1, fcn_sub2 ); 

    auto function1 = fcn_sub_OR1;
    auto function2 = fcn_sub_AND1;
    kitty::dynamic_truth_table func (cut_list.size()+1);
    func._bits[0] = uint64_t(uint64_t(function1._bits[0]) << uint32_t(pow(2,(cut_list.size())))) + function2._bits[0];

    uint32_t cut1_index[5] = {0, 1, 2, 3, 4};
    uint32_t cut2_index[5] = {0, 1, 2, 3, 4};
    bool child_complement[3] = {0};
    uint32_t isize = cut_list.size();
    kitty::dynamic_truth_table function = combine_carry_LUT_function\
      (isize, isize, isize, function1, function2, cut1_index, cut2_index, \
      child_complement);
    
    ntk.set_cell_function(ntk.index_to_node(index), function);

    return true;

  }

  bool set_direct_input (uint32_t cindex, uint32_t index) {
    //std::cout << index << "(direct):";
    std::vector<uint32_t> cut_list;
    carry_cut_list[index].clear();
    for (auto icut: cuts.cuts(index)[0]) {
      if (icut != cindex) {
        cut_list.push_back(icut);
        //std::cout << " " << icut;
        carry_cut_list[index].push_back(icut);
      }
    } 
    carry_cut_list[index].push_back(cindex);

    //std::cout << "\n";
    auto fcn_xsub1 = set_tt_const(index, ntk.get_children(index, 0), cut_list, cindex, 0);
    auto fcn_xsub2 = set_tt_const(index, ntk.get_children(index, 1), cut_list, cindex, 0);
    auto fcn_xsub3 = set_tt_const(index, ntk.get_children(index, 2), cut_list, cindex, 0);
    auto fcn_sub1 = ntk.is_complemented_children(index,0) ? ~fcn_xsub1:fcn_xsub1; 
    auto fcn_sub2 = ntk.is_complemented_children(index,1) ? ~fcn_xsub2:fcn_xsub2;
    auto fcn_sub3 = ntk.is_complemented_children(index,2) ? ~fcn_xsub3:fcn_xsub3;
    auto fcn_sub_OR1 = kitty::binary_or ( fcn_sub1, fcn_sub2 ); 
    auto fcn_sub_OR2 = kitty::binary_or ( fcn_sub1, fcn_sub3 ); 
    auto fcn_sub_OR3 = kitty::binary_or ( fcn_sub2, fcn_sub3 ); 
    auto fcn_sub_AND1 = kitty::binary_and ( fcn_sub1, fcn_sub2 ); 
    auto fcn_sub_AND2 = kitty::binary_and ( fcn_sub1, fcn_sub3 ); 
    auto fcn_sub_AND3 = kitty::binary_and ( fcn_sub2, fcn_sub3 ); 

    auto fcn_top_OR = set_tt_const(index, index, cut_list, cindex, 0);
    auto fcn_top_AND = set_tt_const(index, index, cut_list, cindex, 1);

    auto function1 = fcn_sub_OR1;
    auto function2 = fcn_sub_OR1;
    if ( fcn_sub_OR1 == fcn_top_OR && fcn_sub_AND1 == fcn_top_AND ) {
      function2 = fcn_sub_OR1;
      function1 = fcn_sub_AND1;
    } else if ( fcn_sub_OR1 == fcn_top_AND && fcn_sub_AND1 == fcn_top_OR ) {
      function1 = fcn_sub_OR1;
      function2 = fcn_sub_AND1;
    } else if ( fcn_sub_OR2 == fcn_top_OR && fcn_sub_AND2 == fcn_top_AND ) {
      function2 = fcn_sub_OR2;
      function1 = fcn_sub_AND2;
    } else if ( fcn_sub_OR2 == fcn_top_AND && fcn_sub_AND2 == fcn_top_OR ) {
      function1 = fcn_sub_OR2;
      function2 = fcn_sub_AND2;
    } else if ( fcn_sub_OR3 == fcn_top_OR && fcn_sub_AND3 == fcn_top_AND ) {
      function2 = fcn_sub_OR3;
      function1 = fcn_sub_AND3;
    } else if ( fcn_sub_OR3 == fcn_top_AND && fcn_sub_AND3 == fcn_top_OR ) {
      function1 = fcn_sub_OR3;
      function2 = fcn_sub_AND3;
    } else return false; 

    kitty::dynamic_truth_table func (cut_list.size()+1);
    func._bits[0] = uint64_t(uint64_t(function1._bits[0]) << uint32_t(pow(2,(cut_list.size())))) + function2._bits[0];

    uint32_t cut1_index[5] = {0, 1, 2, 3, 4};
    uint32_t cut2_index[5] = {0, 1, 2, 3, 4};
    bool child_complement[3] = {0};
    uint32_t isize = cut_list.size();
    kitty::dynamic_truth_table function = combine_carry_LUT_function\
      (isize, isize, isize, function1, function2, cut1_index, cut2_index, \
      child_complement);
    
    ntk.set_cell_function(ntk.index_to_node(index), function);

    std::cout << "\t\t" << index << " ";
    kitty::print_hex(func);
    std::cout << ": "; 
    kitty::print_hex(function);
    std::cout << ": ";
    kitty::print_hex (function1);
    std::cout << " "; 
    kitty::print_hex (function2);
    std::cout << "\n";
    cut_list.push_back(cindex);
    auto fcn = flip_carry_user_tt(index,index,cut_list); 
    std::cout << "\t\t" << "original ";
    kitty::print_hex(fcn);
    std::cout << "\n";
    assert(func == fcn);
    return true;

  }

  void characterize_LUTs_helper (uint32_t max_counter, uint32_t index, uint32_t length, uint32_t& path_counter) {

    max_counter++;
    if (max_counter >= MAX_SEARCH || ntk.is_constant(ntk.index_to_node(index))) {
      return;
    } else if (delays[index] == LUT_DELAY) {
      path_counter++;
      return;
    }

    // Find the max delay 
    for ( auto leaf : cuts.cuts(index)[0] ) {
      if (delays[leaf] == delays[index]-LUT_DELAY) {
        characterize_LUTs_helper(max_counter, leaf, length, path_counter);
      }
    }
  }

  void characterize_LUTs (void) {

    uint32_t total_path_counter = 0;
    uint32_t critical_path_counter = 0;

    ntk.foreach_po( [&]( auto const& s ) {
      const auto node = ntk.get_node(s);
      const auto index = ntk.node_to_index(node);
       
      if (ntk.is_pi(node) || ntk.is_constant(node)) return;

      uint32_t path_counter = 0;
      uint32_t max_counter = 0;
      if (delays[index] >= delay-LUT_DELAY) {
        std::cout << "For target index " << index << "(" << delays[index] << "): \n";
        characterize_LUTs_helper (max_counter, index, 1, path_counter);
        critical_path_counter += path_counter;
      }

      total_path_counter += path_counter;

    } );

    std::cout << "This benchmark has " << total_path_counter << " paths from all outputs and "
      << critical_path_counter << " paths from critical outputs.\n";
  }

  void count_mappable_nodes (uint32_t n_crit, uint32_t mig_path_len) {
    uint32_t node_count = 0;
    uint32_t mappable_node_count = 0;
    uint32_t direct_input_mappable = 0;
    ntk.foreach_node( [&]( auto n, auto ) {
      uint32_t index = ntk.node_to_index(n);
      uint32_t num_path = 0;
      uint32_t num_crit_path = 0;
      uint32_t direct_input = 0;
      if ( ntk.is_constant( n ) || ntk.is_pi( n ) ) {
        return;
      }
      //std::cout << index << ":\n";
      for ( auto leaf: cuts.cuts(index)[0] ) {
        if (delays[leaf] == delays[index] - LUT_DELAY) {
          num_crit_path++;
          if ( count_path_to_node (index, index, leaf, 0) == 1  && length_to_node(index, index, leaf) <= mig_path_len ) {
            num_path++;
            //std::cout << "\t" << leaf << "*:" <<count_path_to_node (index, index, leaf, 0) << " "<< length_to_node(index, index, leaf) << "\n";
          }
          if ( check_direct_input (leaf, index) ) {
            direct_input++;
            //std::cout << "\t" << leaf << "^:\n";
          }

          if ( length_to_node(index, index, leaf) > 3 && check_shorter_path(leaf, index, index)) {
            std::cout << index << "->" << leaf << ": " << length_to_node(index, index, leaf) << "\n";
          }
          //if ( count_path_to_node (index, index, leaf, 0) != 1 && check_direct_input(leaf,index) ) {
            //std::cout << index << "->" << leaf << ": " << count_path_to_node (index, index, leaf, 0) << "\n";
            //print_lut(index,index);
            //set_direct_input(leaf,index);
          //}   
          //if ( count_path_to_node (index, index, leaf, 0) == 1 && length_to_node(index, index, leaf) > mig_path_len && check_direct_input(leaf,index) ) {
            //std::cout << index << "->" << leaf << ": " << length_to_node(index, index, leaf) << "\n";
            //print_lut(index,index);
            //set_direct_input(leaf,index);
          //}   
        }
      }
      //std::cout << "\n";
      if (num_crit_path <= n_crit && num_path == num_crit_path) mappable_node_count++;
      if (num_crit_path <= n_crit && direct_input == num_crit_path) direct_input_mappable++;
      if (num_crit_path <= n_crit && num_path != direct_input) {
        //std::cout << index << "\n";
      }
      node_count++;
    });
    std::cout << "There are " << node_count << " nodes and " << mappable_node_count << " nodes are mappable.\n";
    std::cout << "Also, " << direct_input_mappable << " nodes might be mappable\n";
  }
  
  bool is_dependent (uint32_t index, uint32_t cindex, uint32_t curr_index) {
    if (ntk.is_constant(curr_index)) return false;
    //std::cout << "\t\t\t" << curr_index << "\n";
    if (cindex == curr_index) {
      return true;
    }
    for (auto icut: cuts.cuts(index)[0]) {
      if (curr_index == icut) return false;
    }

    for (uint32_t i = 0; i < ntk.fanin_size(curr_index); i++) {
      auto nchild = ntk.get_children(curr_index,i);
      if (ntk.is_constant(nchild) || ntk.is_pi(nchild)) continue;
      if (is_dependent(index, cindex, nchild)) return true;
    }

    return false;
  }

  bool check_shorter_path (uint32_t cindex, uint32_t curr_index, uint32_t index) {
    // check index is not dependent on cindex
    auto child0 = ntk.get_children(curr_index, 0);
    auto child1 = ntk.get_children(curr_index, 1);
    auto child2 = ntk.get_children(curr_index, 2);
    auto indep0 = !is_dependent(index, cindex, child0);
    auto indep1 = !is_dependent(index, cindex, child1);
    auto indep2 = !is_dependent(index, cindex, child2);
    if (indep0 && indep1) {
      //std::cout << "\t1." << child2 << ":\n";
      if (!ntk.is_constant(child2) && child2 != cindex && check_shorter_input (cindex, child2, index)) return true;//std::cout << "HERE\n";
    }
    if (indep0 && indep2) {
      //std::cout << "\t2." << child1 << ":\n";
      if (!ntk.is_constant(child1) && child1 != cindex && check_shorter_input (cindex, child1, index)) return true; //std::cout << "HERE\n";
    }
    if (indep1 && indep2) {
      //std::cout << "\t3." << child0 << ":\n";
      if (!ntk.is_constant(child0) && child0 != cindex && check_shorter_input (cindex, child0, index)) return true; //std::cout << "HERE\n";
    }
    return false;
  }
   bool check_shorter_input (uint32_t cindex, uint32_t index, uint32_t cut_index) {
    //std::cout << index << " " << cindex << ":";
    std::vector<uint32_t> cut_list;
    for (auto icut: cuts.cuts(cut_index)[0]) {
      if (icut != cindex) {
        cut_list.push_back(icut);
        //std::cout << " " << icut;
      }
      if (icut == index) return false;
    } 
    if (cut_list.size() > 5) return false;
    //std::cout << "\n";
    auto fcn_xsub1 = set_tt_const(index, ntk.get_children(index, 0), cut_list, cindex, 0);
    auto fcn_xsub2 = set_tt_const(index, ntk.get_children(index, 1), cut_list, cindex, 0);
    auto fcn_xsub3 = set_tt_const(index, ntk.get_children(index, 2), cut_list, cindex, 0);
    auto fcn_sub1 = ntk.is_complemented_children(index,0) ? ~fcn_xsub1:fcn_xsub1; 
    auto fcn_sub2 = ntk.is_complemented_children(index,1) ? ~fcn_xsub2:fcn_xsub2;
    auto fcn_sub3 = ntk.is_complemented_children(index,2) ? ~fcn_xsub3:fcn_xsub3;
    auto fcn_sub_OR1 = kitty::binary_or ( fcn_sub1, fcn_sub2 ); 
    auto fcn_sub_OR2 = kitty::binary_or ( fcn_sub1, fcn_sub3 ); 
    auto fcn_sub_OR3 = kitty::binary_or ( fcn_sub2, fcn_sub3 ); 
    auto fcn_sub_AND1 = kitty::binary_and ( fcn_sub1, fcn_sub2 ); 
    auto fcn_sub_AND2 = kitty::binary_and ( fcn_sub1, fcn_sub3 ); 
    auto fcn_sub_AND3 = kitty::binary_and ( fcn_sub2, fcn_sub3 ); 

    auto fcn_top_OR = set_tt_const(index, index, cut_list, cindex, 0);
    auto fcn_top_AND = set_tt_const(index, index, cut_list, cindex, 1);

    auto function1 = fcn_sub_OR1;
    auto function2 = fcn_sub_OR1;

    /*std::cout << "\t\t" << index << ":";
    kitty::print_hex(fcn_top_AND);
    std::cout << " "; 
    kitty::print_hex(fcn_top_OR);
    std::cout << ": "; 
    kitty::print_hex(fcn_sub_OR1);
    std::cout << " "; 
    kitty::print_hex(fcn_sub_OR2);
    std::cout << " "; 
    kitty::print_hex(fcn_sub_OR3);
    std::cout << " "; 
    kitty::print_hex(fcn_sub_AND1);
    std::cout << " "; 
    kitty::print_hex(fcn_sub_AND2);
    std::cout << " "; 
    kitty::print_hex(fcn_sub_AND3);
    std::cout << "\n"; 
*/
    if ( fcn_sub_OR1 == fcn_top_OR && fcn_sub_AND1 == fcn_top_AND ) {
      function2 = fcn_sub_OR1;
      function1 = fcn_sub_AND1;
    } else if ( fcn_sub_OR1 == fcn_top_AND && fcn_sub_AND1 == fcn_top_OR ) {
      function1 = fcn_sub_OR1;
      function2 = fcn_sub_AND1;
    } else if ( fcn_sub_OR2 == fcn_top_OR && fcn_sub_AND2 == fcn_top_AND ) {
      function2 = fcn_sub_OR2;
      function1 = fcn_sub_AND2;
    } else if ( fcn_sub_OR2 == fcn_top_AND && fcn_sub_AND2 == fcn_top_OR ) {
      function1 = fcn_sub_OR2;
      function2 = fcn_sub_AND2;
    } else if ( fcn_sub_OR3 == fcn_top_OR && fcn_sub_AND3 == fcn_top_AND ) {
      function2 = fcn_sub_OR3;
      function1 = fcn_sub_AND3;
    } else if ( fcn_sub_OR3 == fcn_top_AND && fcn_sub_AND3 == fcn_top_OR ) {
      function1 = fcn_sub_OR3;
      function2 = fcn_sub_AND3;
    } else return false; 

    kitty::dynamic_truth_table function (cut_list.size()+1);
    function._bits[0] = uint64_t(uint64_t(function1._bits[0]) << uint32_t(pow(2,(cut_list.size())))) + function2._bits[0];
/*
    std::cout << "\t\t" << index << " ";
    kitty::print_hex(function);
    std::cout << ": "; 
    kitty::print_hex (function1);
    std::cout << " "; 
    kitty::print_hex (function2);
    std::cout << "\n";*/
    cut_list.push_back(cindex);
    auto fcn = flip_carry_user_tt(index,index,cut_list); 
    /*std::cout << "\t\t" << "original ";
    kitty::print_hex(fcn);
    std::cout << "\n";*/
    assert(function == fcn);
    return true;

  }
  bool check_direct_input (uint32_t cindex, uint32_t index) {
    //std::cout << index << " " << cindex << ":";
    std::vector<uint32_t> cut_list;
    for (auto icut: cuts.cuts(index)[0]) {
      if (icut != cindex) {
        cut_list.push_back(icut);
        //std::cout << " " << icut;
      }
    } 
    if (cut_list.size() > 5) return false;
    //std::cout << "\n";
    auto fcn_xsub1 = set_tt_const(index, ntk.get_children(index, 0), cut_list, cindex, 0);
    auto fcn_xsub2 = set_tt_const(index, ntk.get_children(index, 1), cut_list, cindex, 0);
    auto fcn_xsub3 = set_tt_const(index, ntk.get_children(index, 2), cut_list, cindex, 0);
    auto fcn_sub1 = ntk.is_complemented_children(index,0) ? ~fcn_xsub1:fcn_xsub1; 
    auto fcn_sub2 = ntk.is_complemented_children(index,1) ? ~fcn_xsub2:fcn_xsub2;
    auto fcn_sub3 = ntk.is_complemented_children(index,2) ? ~fcn_xsub3:fcn_xsub3;
    auto fcn_sub_OR1 = kitty::binary_or ( fcn_sub1, fcn_sub2 ); 
    auto fcn_sub_OR2 = kitty::binary_or ( fcn_sub1, fcn_sub3 ); 
    auto fcn_sub_OR3 = kitty::binary_or ( fcn_sub2, fcn_sub3 ); 
    auto fcn_sub_AND1 = kitty::binary_and ( fcn_sub1, fcn_sub2 ); 
    auto fcn_sub_AND2 = kitty::binary_and ( fcn_sub1, fcn_sub3 ); 
    auto fcn_sub_AND3 = kitty::binary_and ( fcn_sub2, fcn_sub3 ); 

    auto fcn_top_OR = set_tt_const(index, index, cut_list, cindex, 0);
    auto fcn_top_AND = set_tt_const(index, index, cut_list, cindex, 1);

    auto function1 = fcn_sub_OR1;
    auto function2 = fcn_sub_OR1;

    /*std::cout << "\t\t" << index << ":";
    kitty::print_hex(fcn_top_AND);
    std::cout << " "; 
    kitty::print_hex(fcn_top_OR);
    std::cout << ": "; 
    kitty::print_hex(fcn_sub_OR1);
    std::cout << " "; 
    kitty::print_hex(fcn_sub_OR2);
    std::cout << " "; 
    kitty::print_hex(fcn_sub_OR3);
    std::cout << " "; 
    kitty::print_hex(fcn_sub_AND1);
    std::cout << " "; 
    kitty::print_hex(fcn_sub_AND2);
    std::cout << " "; 
    kitty::print_hex(fcn_sub_AND3);
    std::cout << "\n"; 
*/
    if ( fcn_sub_OR1 == fcn_top_OR && fcn_sub_AND1 == fcn_top_AND ) {
      function2 = fcn_sub_OR1;
      function1 = fcn_sub_AND1;
    } else if ( fcn_sub_OR1 == fcn_top_AND && fcn_sub_AND1 == fcn_top_OR ) {
      function1 = fcn_sub_OR1;
      function2 = fcn_sub_AND1;
    } else if ( fcn_sub_OR2 == fcn_top_OR && fcn_sub_AND2 == fcn_top_AND ) {
      function2 = fcn_sub_OR2;
      function1 = fcn_sub_AND2;
    } else if ( fcn_sub_OR2 == fcn_top_AND && fcn_sub_AND2 == fcn_top_OR ) {
      function1 = fcn_sub_OR2;
      function2 = fcn_sub_AND2;
    } else if ( fcn_sub_OR3 == fcn_top_OR && fcn_sub_AND3 == fcn_top_AND ) {
      function2 = fcn_sub_OR3;
      function1 = fcn_sub_AND3;
    } else if ( fcn_sub_OR3 == fcn_top_AND && fcn_sub_AND3 == fcn_top_OR ) {
      function1 = fcn_sub_OR3;
      function2 = fcn_sub_AND3;
    } else return false; 

    kitty::dynamic_truth_table function (cut_list.size()+1);
    function._bits[0] = uint64_t(uint64_t(function1._bits[0]) << uint32_t(pow(2,(cut_list.size())))) + function2._bits[0];

    /*std::cout << "\t\t" << index << " ";
    kitty::print_hex(function);
    std::cout << ": "; 
    kitty::print_hex (function1);
    std::cout << " "; 
    kitty::print_hex (function2);
    std::cout << "\n";*/
    cut_list.push_back(cindex);
    auto fcn = flip_carry_user_tt(index,index,cut_list); 
    /*std::cout << "\t\t" << "original ";
    kitty::print_hex(fcn);
    std::cout << "\n";*/
    assert(function == fcn);
    return true;

  }

  kitty::dynamic_truth_table set_tt_const ( uint32_t index, uint32_t curr_index, auto cut_leaf_list, uint32_t const_leaf, bool val ) {
 
    kitty::dynamic_truth_table tt_new (cut_leaf_list.size());

    if ( ntk.is_constant(curr_index) ) return tt_new;
    if ( curr_index == const_leaf ) {
      if (val) return ~tt_new;//._bits[0]=0xFFFFFFFFFFFFFFFF;
      else return tt_new;
    }

    uint32_t place = 0;
    for (auto leaf: cut_leaf_list) {
      if (leaf == curr_index ) {
        for (uint64_t k = 0; k < pow(2,cut_leaf_list.size()); k++) {
          tt_new._bits[0] += uint64_t(get_bit(k, place) << k);
        }
        return tt_new;
      }
      place++;
    } 
 
    auto n = ntk.index_to_node(curr_index);

    uint32_t j = 0; 
    std::vector<kitty::dynamic_truth_table> tt_child(3);
    ntk.foreach_fanin( n, [&]( auto const& f ) {
      auto n_child = ntk.get_node(f);
      tt_child[j] = set_tt_const(index, n_child, cut_leaf_list, const_leaf, val);
      j++;  
    });
  
    tt_new = ntk.compute(n, tt_child.begin(), tt_child.end());
    return tt_new;
    
  }

  // How many critical POs (not paths) does the node connect to eventually
  void set_critical_node_helper (uint32_t index, uint32_t& max_counter) {

    max_counter++;
    if (max_counter >= MAX_SEARCH || ntk.is_pi(ntk.index_to_node(index)) || ntk.is_constant(ntk.index_to_node(index))) return;

    num_critical_path_through_node[index]++;

    // Find the max delay 
    for ( auto leaf : cuts.cuts(index)[0] ) {
      if (delays[leaf] == delays[index]-LUT_DELAY) {
        set_critical_node_helper(leaf, max_counter);
      }
    }
  }

  void set_critical_node (void) {

    ntk.foreach_po( [&]( auto const& s ) {
      const auto node = ntk.get_node(s);
      const auto index = ntk.node_to_index(node);
       
      if (ntk.is_pi(node) || ntk.is_constant(node)) return;

      uint32_t max_counter = 0; 
      if (delays[index] >= delay) {
        set_critical_node_helper(index, max_counter);
        num_critical_map_refs[index]++;
      }
    } );

    for ( auto it = top_order.rbegin(); it != top_order.rend(); ++it )  {
      /* skip constants and PIs (TODO: stop earlier) */
      if ( ntk.is_constant( *it ) || ntk.is_pi( *it ) )
        continue;

      const auto index = ntk.node_to_index( *it );
      if ( num_critical_map_refs[index] == 0 )
        continue;

      for ( auto leaf : cuts.cuts( index )[0] ) {
        if (delays[leaf] == delays[index]-LUT_DELAY)
          num_critical_map_refs[leaf]++;
      }
    }

    /*std::cout << "Critical nodes are:\n";
    for (int32_t i = num_critical_path_through_node.size() - 1; i >= 0; i--) {
      if (num_critical_path_through_node[i] > 0)
        std::cout << i << "(" << delays[i] << ") = " << num_critical_map_refs[i] << " " << map_refs[i] << " " << num_critical_path_through_node[i] << "\n"; 
    }*/
  }
 
  bool path_selection (std::vector<node<Ntk>>& path_for_carry_chain, uint32_t delay_offset) {

    for (uint32_t curr_offset = 0; curr_offset <= delay_offset; curr_offset+=LUT_DELAY) {
      for ( auto it = top_order.rbegin(); it != top_order.rend(); ++it ) {
      //ntk.foreach_po( [&]( auto const& s ) {
        const auto node = *it;
        //const auto node = ntk.get_node(s);
        const auto index = ntk.node_to_index(node);

        if (!path_for_carry_chain.empty()) break;
        // Node with worst delay
        if (delays[index] >= delay-curr_offset && delay != 0 && map_refs[index] > 0  && !ntk.is_pi(node)) {
          std::cout << "Found target index " << index << "(" << delays[index] << ")\n";

          // Try to place a path starting from this node
          if (find_deepest_LUT(path_for_carry_chain,index, 0)) {
            if (!is_a_carry_node(index)) {
              std::cout << "adding " << index << "\n";
              path_for_carry_chain.push_back(index);
              //carry_nodes[index] += 1;
            }
            // Only place one path at a time
            //return;
            break;
          } else {
            path_for_carry_chain.clear();
          }
        }
      }// );
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
              // ALLOW THIS IF CARRY IN WILL BE CONNECTED THROUGH LUT
              if (carry_i == 1) std::cout << n << " IS BAD CARRY IN\n"; 
              //assert(0);
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
          //carry_truth_table (n, carry_driver_nodes[index]);
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

    // Traverse through the list backwards to avoid carry in with dependencies on path after
    // This can happen due to multiple paths being placed avoiding mapped carry in
    for (int i = carry_paths.size()-1; i >= 0; i--) {
      auto carry_path  = carry_paths[i];
 
      for (uint32_t carry_i = 0; carry_i < carry_path.size(); carry_i++) {
        
        auto carry_node = carry_path[carry_i];
        auto carry_child_node = carry_driver_nodes[carry_node];
        
        // Only consider flipping carry in signal if it has a driver (from another carry chain)
        if (carry_child_node == 0) continue;

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
          }
          if (!is_a_carry_node(child_node) && mapped_to_5LUTb[carry_node] != 0 && mapped_to_5LUTb[carry_node] == child_node) {
            mapped_to_5LUT_complemented[child_node] = ntk.is_complemented_children(carry_node,i);
          }
        }
        if (carry_i == carry_path.size()-1) {
          ntk.foreach_po( [this]( auto s ) {
            const auto index = ntk.node_to_index( ntk.get_node( s ) );
            if (is_a_carry_node(ntk.get_node(s)) && ntk.is_complemented(s)) {
              carry_nodes[ntk.index_to_node(index)] = 0;
              
              std::cout << "Undoing carry for the last PO node if inverted\n";

              //delays[index] += LUT_DELAY;
              //std::cout << "Updating this particular delay " << index << "(" << delays[index] << ")\n";
            }
            delay = std::max( delay, delays[index] );
          });

        }
      }
    }

    delay = 0;
    ntk.foreach_po( [this]( auto s ) {
      const auto index = ntk.node_to_index( ntk.get_node( s ) );
      if (is_a_carry_node(ntk.get_node(s)) && ntk.is_complemented(s)) {
        carry_nodes[ntk.index_to_node(index)] = 0;
        
        //std::cout << "Undoing carry for the last PO node if inverted\n";

        delays[index] += LUT_DELAY;
        std::cout << "Updating this particular delay " << index << "(" << delays[index] << ")\n";
      }
      delay = std::max( delay, delays[index] );
    });
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
      //std::cout << index << ":";
      if ( ntk.is_constant( n ) || ntk.is_pi( n ) ) {
        /* all terminals have flow 1.0 */
        flow_refs[index] = 1.0f;
      } else {
        flow_refs[index] = static_cast<float>( ntk.fanout_size( n ) );
      }
      flows[index] = cuts.cuts( index )[0]->data.flow;
      delays[index] = cuts.cuts( index )[0]->data.delay*LUT_DELAY;
      
    } );

    // Carry chain related initialization
    tmp_size = ntk.size();
    for (uint32_t i = 0; i < ntk.size(); i++) {
      remapped_nodes[i].resize(tmp_size, -1);
    }
  }

  ///////////////////////////////////////////////////////////////
  // Map to LUTs before carry nodes 
  ///////////////////////////////////////////////////////////////

  void xilinx_compute_carry_mapping (std::vector<node<Ntk>> path_for_carry_chain)
  {

    node<Ntk> n_carryin, n_first;

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

        cost = cost_of_5LUT_cuts(i_child[0], i_child[1], cut_i_1, cut_i_2);

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
    //bool match = false;

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

        //bool newcheck = check_cut_legality(i1_child[0], i1_child[1], i2_child[0], i2_child[1],
        //  cut_i1_1, cut_i1_2, cut_i2_1, cut_i2_2);

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
        double cost = cost_of_cuts(i1_child[0], i1_child[1],
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
          //print_cut_cost (i1, i1_child[0], i1_child[1], cut_i1_1, cut_i1_2, cost); 
          //print_cut_cost (i2, i2_child[0], i2_child[1], cut_i2_1, cut_i2_2, cost); 
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
      uint32_t new_delay = 0; 
      //std::cout << index << ":";
      if (ntk.is_pi(n) || ntk.is_constant(n)) {
        //std::cout << "PI or CONST\n";
        continue;
      } else if (is_a_carry_node(n)) {
        //new_delay = delays[carry_driver_nodes[index]] + CARRY_DELAY;
        //std::cout << carry_driver_nodes[index] << "(" << new_delay << "),";
        for (auto leaf: carry_cut_list[index]) {
          uint32_t leaf_delay = delays[leaf];
        
          // leafs from carry node get lut + adder delay
          if (carry_driver_nodes[index] == leaf) leaf_delay += CARRY_DELAY;
          else  leaf_delay += LUT_ADDER_DELAY; 
          //std::cout << leaf << "(" << leaf_delay << "),";
          if (leaf_delay > new_delay) {
            new_delay = leaf_delay;
          }
        }
      } else { 
        
        for (auto leaf: cuts.cuts(index)[0]) {
          uint32_t leaf_delay = delays[leaf];
          leaf_delay += LUT_DELAY;
          // leafs from carry node get lut + adder delay
          //if (carry_driver_nodes[index] == leaf) leaf_delay += CARRY_DELAY;
          //else  leaf_delay += LUT_ADDER_DELAY; 
          //std::cout << leaf << "(" << leaf_delay << "),";
          if (leaf_delay > new_delay) {
            new_delay = leaf_delay;
          }
        }
        //std::cout << "(" << new_delay << "),";
      }
      delays[index] = new_delay;
      //std::cout << "\n";
    } 

    delay = 0;
    ntk.foreach_po( [this]( auto s ) {
      const auto index = ntk.node_to_index( ntk.get_node( s ) );
      delay = std::max( delay, delays[index] );
    });
  }

  void update_area() {

    std::vector<uint32_t> temp_map_refs(ntk.size(), 0);
    /* compute current delay */
    ntk.foreach_po( [&]( auto s ) {
      const auto index = ntk.node_to_index( ntk.get_node( s ) );
      temp_map_refs[index]++;
    });

    /* compute current area and update mapping refs */
    area = 0;
    for ( auto it = top_order.rbegin(); it != top_order.rend(); ++it )
    {
      /* skip constants and PIs (TODO: stop earlier) */
      if ( ntk.is_constant( *it ) || ntk.is_pi( *it ) )
        continue;

      const auto index = ntk.node_to_index( *it );
      if ( temp_map_refs[index] == 0 )
        continue;
     
      if (is_a_carry_node( *it )) {
        for ( auto leaf : carry_cut_list[index] )
        {
          temp_map_refs[leaf]++;
        }
  
      } else {
        for ( auto leaf : cuts.cuts( index )[0] )
        {
          temp_map_refs[leaf]++;
        }
      }
      area++;
    }
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
      std::cout << "set " << index << "\n";
    });

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
     
      std::cout << "because of " << index << ": ";
 
      if (is_a_carry_node( *it )) {
        for ( auto leaf : carry_cut_list[index] )
        {
          map_refs[leaf]++;
          std::cout << leaf << " ";
        }
  
      } else {
        for ( auto leaf : cuts.cuts( index )[0] )
        {
          map_refs[leaf]++;
          std::cout << leaf << " ";
        }
      }
      std::cout << "\n";
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

 
    if ( (( ntotal_index1 > 3 ) && ((ntotal_index2 - nmatches) > 2)) || \
         (( ntotal_index2 > 3 ) && ((ntotal_index1 - nmatches) > 2)))
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
    assert(curr_LUT <= 2);
  }

  void get_children_flipped (node<Ntk>child_index[2], uint32_t index, uint32_t carryin) {

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
    assert(curr_LUT <= 2);
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

  kitty::dynamic_truth_table combine_carry_LUT_function (uint32_t n_total, uint32_t n_cut1, uint32_t n_cut2,
      kitty::dynamic_truth_table cut1_function, kitty::dynamic_truth_table cut2_function, 
      uint32_t cut1_index[5], uint32_t cut2_index[5], bool child_complement[3]) {

    // Create a truth table
    uint32_t lut_input_size = n_total + 1; 
    kitty::dynamic_truth_table function (lut_input_size);
    function._bits[0]= 0x0000000000000000;

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
        for (uint32_t j = 0; j < n_cut1; j++) {
          cut1_i += (get_bit(i, cut1_index[j]) << j);
        } 
        uint32_t cut2_i = 0;
        for (uint32_t j = 0; j < n_cut2; j++) {
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
        for (uint32_t j = 0; j < n_cut1; j++) {
          cut1_i += (get_bit(i, cut1_index[j]) << j);
        } 
        uint32_t cut2_i = 0;
        for (uint32_t j = 0; j < n_cut2; j++) {
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
    return function;
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

    auto cut1_function_bu = cuts.truth_table( cut1 ); 
    auto cut2_function_bu = cuts.truth_table( cut2 ); 
    
    auto cut1_function = flip_carry_user_tt (cindex_1, cindex_1, cuts.cuts(cindex_1)[a] );
    auto cut2_function = flip_carry_user_tt (cindex_2, cindex_2, cuts.cuts(cindex_2)[b] );

    return combine_carry_LUT_function (n_total, cut1.size(), cut2.size(), cut1_function, cut2_function, 
      cut1_index, cut2_index, child_complement);
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
        std::cout << "Error: " << n << " does not have child " << ntk.get_node(f) << "\n";
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

    //auto n = ntk.index_to_node(index);

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
    //uint32_t index_child_cuts[2][5] = {0}; 
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
  double cost_of_5LUT_cuts(uint32_t index_1, uint32_t index_2, int32_t i1, uint32_t i2) {

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

  double cost_of_cuts( uint32_t index1_1, uint32_t index1_2,
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
      cost *= (1 + nshared1/ntotal1 + nshared2/ntotal2 + nshared/(ntotal1+ntotal2));
      
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
    //uint32_t nmatch_1_1_2_1 = number_of_match_in_cuts(cut_1_1, cut_2_1);
    //uint32_t nmatch_1_1_2_2 = number_of_match_in_cuts(cut_1_1, cut_2_2);
    //uint32_t nmatch_1_2_2_1 = number_of_match_in_cuts(cut_1_2, cut_2_1);
    //uint32_t nmatch_1_2_2_2 = number_of_match_in_cuts(cut_1_2, cut_2_2);
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
  void compute_mapping( bool delay, bool baseline )
  {
    for ( auto const& n : top_order )
    {
      if ( ntk.is_constant( n ) || ntk.is_pi( n ) || is_a_carry_node( n ) )
        continue;
      compute_best_cut<ELA>( ntk.node_to_index( n ), delay, baseline );
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
    uint32_t mappable_area = 0;
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

    return {(flow + cut_area( cut )), time+LUT_DELAY};
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
  void compute_best_cut( uint32_t index, bool delay, bool baseline )
  {
    constexpr auto mf_eps{0.005f};

    float flow;
    uint32_t time{0};
    int32_t best_cut{-1};
    float best_flow{std::numeric_limits<float>::max()};
    uint32_t best_time{std::numeric_limits<uint32_t>::max()};
    int32_t cut_index{-1};
    uint32_t cost = 0;
    uint32_t best_cost = 0;

    //std::cout << "index " << index << "\n";
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

      uint32_t path_factor = 0;
      uint32_t n_crit = 0;
      for ( auto leaf : *cut ) {
        if ( (time - LUT_DELAY) == delays[leaf] ) {
          n_crit++;
          //if (count_path_to_node (index,index,leaf,cut_index) == 1) path_factor++;
        } 
      } 
     
      if (baseline) cost = 0;
      //else if (n_crit == 1) {
      //  cost = 0;
        //if (delay) flow *= 0.9;
        //else time -= LUT_DELAY;
      //} else if (n_crit == 2) cost = 1;
      //else cost = 2;
      //else cost = (n_crit - path_factor);
      /*auto temp_flow = 0; 
      if (cost != 0) {
        temp_flow = flow + 1.0;
      }*/

      if ( delay ) {
        if ( best_cut == -1 || best_time > time || ( best_time == time &&  best_flow > flow + mf_eps ) 
        || ( best_time == time && best_flow > flow - mf_eps && best_cost > cost ) )
        {
          best_cut = cut_index;
          best_flow = flow;
          best_time = time;
          best_cost = cost;
        }
      } else {
        if ( best_cut == -1 || best_flow > flow + mf_eps || ( best_flow > flow - mf_eps && best_time > time ) 
        || ( best_flow > flow - mf_eps && best_time == time && best_cost > cost ) )
        {
          best_cut = cut_index;
          best_flow = flow;
          best_time = time;
          best_cost = cost;
          //std::cout << "*";
        }
      }
      //std::cout << index << " " << cut_index << ": " << n_crit << " " << path_factor << " " << flow << " " << time << " -> " << cost << "\n";
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

        
        const auto index = ntk.node_to_index( n );

        auto updated_n = n;
        auto updated_index = index;
        if (index < tmp_size && remapped_nodes[index][index] > -1) {
          updated_n = ntk.get_node(tmp_signals[remapped_nodes[index][index]]);
          updated_index = ntk.node_to_index(n);
        }

        if (ps.verbose && ps.verbosity > 3) std::cout << "Node* " << updated_n;

        // Add node as carry
        std::vector<node<Ntk>> nodes;
        ntk.add_to_carry_mapping ( updated_n, ntk.node_to_index(n_carry) ); 
        
        // Add LUT mapped before carry  
        uint32_t index_c1 = 0;
        uint32_t index_c2 = 0;
        if (!is_a_carry_node(mapped_to_5LUTa[updated_index]) && mapped_to_5LUTa[updated_index] > 0) 
          index_c1 = mapped_to_5LUTa[updated_index];
        if (!is_a_carry_node(mapped_to_5LUTb[updated_index]) && mapped_to_5LUTb[updated_index] > 0)
          index_c2 = mapped_to_5LUTb[updated_index];

        ntk.add_to_carry_LUT_mapping( updated_n, index_c1, index_c2, \
            mapped_to_5LUT_complemented[index_c1], mapped_to_5LUT_complemented[index_c2] );
        if (ps.verbose && ps.verbosity > 3) std::cout << "(" << index_c1 << "/" << index_c2 << ")" << ":";

        for ( auto c: carry_cut_list[updated_index]) {
          if (c < tmp_size && remapped_nodes[c][c] > -1) 
            c = ntk.get_node(tmp_signals[remapped_nodes[c][c]]);
          nodes.push_back(c);
          if (ps.verbose && ps.verbosity > 3) std::cout << c << " ";
        }
        if (ps.verbose && ps.verbosity > 3) std::cout << "\n";
        if (ps.verbose && ps.verbosity > 3) std::cout << "\t";
        if (ps.verbose && ps.verbosity > 3) kitty::print_hex(ntk.cell_function(n));
        if (ps.verbose && ps.verbosity > 3) std::cout << "\n";
        ntk.add_to_mapping( updated_n, nodes.begin(), nodes.end() );
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
      for ( auto c : cuts.cuts( index ).best() ) {
        if (c < tmp_size && remapped_nodes[c][c] > -1) 
          c = ntk.get_node(tmp_signals[remapped_nodes[c][c]]);
        nodes.push_back(c);
        if (ps.verbose && ps.verbosity > 3) std::cout << c << " ";
      }
      if (ps.verbose && ps.verbosity > 3) std::cout << "\n";
      ntk.add_to_mapping( n, nodes.begin(), nodes.end() );

      // This is translated into LUT functionality
      if constexpr ( StoreFunction )
      {
        if ( !is_a_carry_node(n) ) {
          auto tt_updated = flip_carry_user_tt ( index, index, cuts.cuts(index)[0] );

          ntk.set_cell_function( n, tt_updated );

          if (ps.verbose && ps.verbosity > 3) std::cout << "\t";
          if (ps.verbose && ps.verbosity > 3) kitty::print_hex(ntk.cell_function(n));
          if (ps.verbose && ps.verbosity > 3) std::cout << "\n";
        }
      }
    }

    for (uint32_t i = 0; i < tmp_size; i++) {
      if (remapped_nodes[i][i] > -1) {
        ntk.substitute_node(i, tmp_signals[remapped_nodes[i][i]]);
        std::cout << "substituting " << i << " with " << ntk.get_node(tmp_signals[remapped_nodes[i][i]]) << "\n";
      }
    }
  }

  /*kitty::dynamic_truth_table flip_carry_user_tt ( uint32_t index, uint32_t curr_index, auto cut_leaf_list ) {
 
    kitty::dynamic_truth_table tt_new (cut_leaf_list.size());

    if ( ntk.is_constant(curr_index)) return tt_new;
    //if ( ntk.is_pi(curr_index)) {
    //  for (uint64_t k = 0; k < pow(2,cut_leaf_list.size()); k++) {
    //    tt_new._bits[0] += uint64_t(get_bit(k,0) << k);
    //  }
    //  return tt_new;
    //}

    uint32_t place = 0;
    bool here = false;
    bool flip = false;
    for (auto leaf: cut_leaf_list) {
      if (leaf == curr_index ) {
        if(ps.verbose && ps.verbosity > 4) std::cout << "\tfound leaf " << leaf << "\n";
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
      for (uint64_t k = 0; k < pow(2,cut_leaf_list.size()); k++) {
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
      tt_child[j] = flip_carry_user_tt(index, n_child, cut_leaf_list);
      j++;  
    });
  
    tt_new = ntk.compute(n, tt_child.begin(), tt_child.end());
    return tt_new;
    
  }*/
  kitty::dynamic_truth_table flip_carry_user_tt ( uint32_t index, uint32_t curr_index, auto cut_leaf_list, bool use_remap = false, bool input_flip = false ) {
    kitty::dynamic_truth_table tt_new (cut_leaf_list.size());

    bool flipped = false;
    if(use_remap && index < tmp_size && curr_index < tmp_size && remapped_nodes[index][curr_index] > -1 ) {
      //std::cout << "from " << curr_index << "(" << remapped_nodes[index][curr_index] << ")";
      curr_index = ntk.node_to_index(ntk.get_node(tmp_signals[remapped_nodes[index][curr_index]]));
      //std::cout << " to " << curr_index;
      flipped = input_flip;
      //std::cout << "\n";
    }

    if ( ntk.is_constant(curr_index) && flipped) return ~tt_new;
    if ( ntk.is_constant(curr_index)) return tt_new;
    /*if ( ntk.is_pi(curr_index)) {
      for (uint64_t k = 0; k < pow(2,cut_leaf_list.size()); k++) {
        tt_new._bits[0] += uint64_t(get_bit(k,0) << k);
      }
      return tt_new;
    }*/

    uint32_t place = 0;
    bool here = false;
    bool flip = false;
    for (auto leaf: cut_leaf_list) {
      if (leaf == curr_index ) {
        if(ps.verbose && ps.verbosity > 4) std::cout << "\tfound leaf " << leaf << "\n";
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
      for (uint64_t k = 0; k < pow(2,cut_leaf_list.size()); k++) {
        //if (flip == true )
        //  tt_new._bits[0] += uint64_t(get_not_bit(k, place) << k);
        //else 
        tt_new._bits[0] += uint64_t(get_bit(k, place) << k);
      }
      if (flip || flipped) {
        return ~tt_new;
      }
      return tt_new;
    }

    auto n = ntk.index_to_node(curr_index);

    uint32_t j = 0; 
    std::vector<kitty::dynamic_truth_table> tt_child(3);
    ntk.foreach_fanin( n, [&]( auto const& f ) {
      auto n_child = ntk.get_node(f);
      tt_child[j] = flip_carry_user_tt(index, n_child, cut_leaf_list, use_remap, input_flip);
      //std::cout << "\t\t" << n_child << ": "; kitty::print_hex(tt_child[j]); std::cout <<"\n";
      j++;  
    });
  
    tt_new = ntk.compute(n, tt_child.begin(), tt_child.end());
    //std::cout << "\t\t" << n << ": "; kitty::print_hex(tt_new); std::cout <<"\n";

    if (flipped) return ~tt_new;
    return tt_new;
    
  }


  ///////////////////////////////////////////////////////////////
  // Rewrite graph for better mapping 
  ///////////////////////////////////////////////////////////////
  bool relevance ( node<Ntk> const& n, uint32_t leaf ) {

    // Get order of children (depends on rand)
    const auto ocs = ordered_children( n, 0 );

    auto index = ntk.node_to_index(n);
    auto curr_tt = flip_carry_user_tt(n,n, cuts.cuts(index)[0]);
    auto new_tt = curr_tt;

    bool mappable = false;

    // Relevance
    if ( auto cand = relevance_candidate ( ocs[0], ocs[1], ocs[2], n, leaf ); cand ) { 
      const auto& [x,y,z] = *cand;
      //print_lut(n,n);
      std::cout << "\t\t\t\tapply relevance to " << n << "\n";

      //auto opt = ntk.create_maj( x, y, z );
      //auto new_node = ntk.get_node(opt);
      //tmp_signals.push_back(opt);
      //remapped_nodes[n][n] = tmp_signals.size()-1; 
      remapped_nodes[n][n] = -2; 
      
      tmp_signals.push_back(y);
      remapped_nodes[n][ntk.get_node(x)] = tmp_signals.size()-1; 
      //print_lut(n,n, true);
      if ( count_path_to_node (n, ntk.get_node(z), leaf, 0, true) == 0) {
        mappable = true;
 
        bool flipped = !(ntk.is_complemented(y) ^ ntk.is_complemented(x));
        std::vector<kitty::dynamic_truth_table> tt_child(3);

        remapped_nodes_flip[n] = flipped;
        if (z == ocs[0]) {
          tt_child[0] = flip_carry_user_tt(n, ntk.get_node(z), cuts.cuts(index)[0], true, flipped);
          tt_child[1] = flip_carry_user_tt(n, ntk.get_node(ocs[1]), cuts.cuts(index)[0]);
          tt_child[2] = flip_carry_user_tt(n, ntk.get_node(ocs[2]), cuts.cuts(index)[0]);
        } else if (z == ocs[1]) {
          tt_child[0] = flip_carry_user_tt(n, ntk.get_node(ocs[0]), cuts.cuts(index)[0]);
          tt_child[1] = flip_carry_user_tt(n, ntk.get_node(z), cuts.cuts(index)[0], true, flipped);
          tt_child[2] = flip_carry_user_tt(n, ntk.get_node(ocs[2]), cuts.cuts(index)[0]);
        } else {
          tt_child[0] = flip_carry_user_tt(n, ntk.get_node(ocs[0]), cuts.cuts(index)[0]);
          tt_child[1] = flip_carry_user_tt(n, ntk.get_node(ocs[1]), cuts.cuts(index)[0]);
          tt_child[2] = flip_carry_user_tt(n, ntk.get_node(z), cuts.cuts(index)[0], true, flipped);
        }
        //std::cout << "\t\tafter: "; kitty::print_hex(tt_child[0]); std::cout <<"\n";
        //std::cout << "\t\tafter: "; kitty::print_hex(tt_child[1]); std::cout <<"\n";
        //std::cout << "\t\tafter: "; kitty::print_hex(tt_child[2]); std::cout <<"\n";
        new_tt = ntk.compute(n, tt_child.begin(), tt_child.end());

      } else assert ( 0 && "Relevance did not remove all paths!");
    } 

    if(curr_tt != new_tt) {
      std::cout << "\t\t\t\t\tbefore: "; kitty::print_hex(curr_tt); std::cout <<"\n";
      std::cout << "\t\t\t\t\tafter: "; kitty::print_hex(new_tt); std::cout <<"\n";
      assert(0 && "After mapping, wrong truth table.\n");
    }

    return mappable;
  }

  using candidate3_t = std::tuple<signal<Ntk>, signal<Ntk>, signal<Ntk>>;
  std::optional<candidate3_t> relevance_candidate ( signal<Ntk> const& a, signal<Ntk> const& b, signal<Ntk> const& c, node<Ntk> n, uint32_t leaf ) 
  {
    // leaf exists in 2 children?
    auto index = ntk.node_to_index(n);
    uint32_t path_in_a = count_path_to_node (index, ntk.get_node(a), leaf, 0, true);
    uint32_t path_in_b = count_path_to_node (index, ntk.get_node(b), leaf, 0, true);
    uint32_t path_in_c = count_path_to_node (index, ntk.get_node(c), leaf, 0, true);
    bool leaf_is_a = ntk.get_node(a) == leaf;
    bool leaf_is_b = ntk.get_node(b) == leaf;
    bool leaf_is_c = ntk.get_node(c) == leaf;

    if (leaf_is_a && ((path_in_b == 0 && path_in_c > 0) || (path_in_c == 0 && path_in_b > 0))){
      //std::cout << "abc is " << path_in_a << path_in_b << path_in_c << "\n";
      return path_in_c > 0 ? candidate3_t{a, b, c} : candidate3_t{a, c, b};
    }
    if (leaf_is_b && ((path_in_a == 0 && path_in_c > 0) || (path_in_c == 0 && path_in_a > 0))){
      //std::cout << "abc is " << path_in_a << path_in_b << path_in_c << "\n";
      return path_in_c > 0 ? candidate3_t{b, a, c} : candidate3_t{b, c, a};
    }
    if (leaf_is_c && ((path_in_b == 0 && path_in_a > 0) || (path_in_a == 0 && path_in_b > 0))){
      //std::cout << "abc is " << path_in_a << path_in_b << path_in_c << "\n";
      return path_in_a > 0 ? candidate3_t{c, b, a} : candidate3_t{c, a, b};
    }
    return std::nullopt;
  }

  bool associativity ( node<Ntk> const& n, uint32_t leaf, uint32_t rand ) {

    // Get order of children (depends on rand)
    const auto ocs = ordered_children( n, rand );

    // Check that the 3rd input has children to consider
    if ( !ntk.is_maj( ntk.get_node( ocs[2] ) ) )
      return false;

    auto index = ntk.node_to_index(n);
    for ( auto cut_index: cuts.cuts(index)[0] ) {
      if (ntk.get_node(ocs[2]) == cut_index) return false;
    } 

    // Get children of 3rd input
    auto ocs2 = ordered_children( ntk.get_node( ocs[2] ), 0 );

    // Propagate inverter
    if ( ntk.is_complemented( ocs[2] ) ) {
      ocs2[0] = !ocs2[0];
      ocs2[1] = !ocs2[1];
      ocs2[2] = !ocs2[2];
    }

    auto curr_tt = flip_carry_user_tt(n,n, cuts.cuts(index)[0]);
    auto new_tt = curr_tt;
    bool mappable = false;

    if ( auto cand = associativity_candidate( ocs[0], ocs[1], ocs2[0], ocs2[1], ocs2[2] ); cand ) {
      const auto& [x, y, z, u, assoc] = *cand;
      auto opt = ntk.create_maj( z, assoc ? u : x, ntk.create_maj( x, y, u ) );
     
      auto new_node = ntk.get_node(opt);
      tmp_signals.push_back(opt);
      remapped_nodes[n][n] = tmp_signals.size()-1; 

      if ( count_path_to_node (n, n, leaf, 0, true) == 1) {
        //std::cout << "\t\t\t\tassociativity 1 " << n << " with " << new_node << "\n";
        mappable = true;
        new_tt = flip_carry_user_tt(n, n, cuts.cuts(index)[0]);
        //print_lut(n,n, false);
        //print_lut(n,n, true);
      }
      if (!mappable)  {
        tmp_signals.pop_back();
        remapped_nodes[n][n] = -1;
      }
    } 

    if ( auto cand = associativity_candidate( ocs[0], ocs[1], ocs2[0], ocs2[2], ocs2[1] ); !mappable && cand ) {
      const auto& [x, y, z, u, assoc] = *cand;
      auto opt = ntk.create_maj( z, assoc ? u : x, ntk.create_maj( x, y, u ) );
     
      auto new_node = ntk.get_node(opt);
      tmp_signals.push_back(opt);
      remapped_nodes[n][n] = tmp_signals.size()-1; 

      if ( count_path_to_node (n, n, leaf, 0, true) == 1) {
        //std::cout << "\t\t\t\tassociativity 2 " << n << " with " << new_node << "\n";
        mappable = true;
        new_tt = flip_carry_user_tt(n, n, cuts.cuts(index)[0]);
        //ntk.set_cell_function(n, new_tt);
        //print_lut(n,n, false);
        //print_lut(n,n, true);
      }
      if (!mappable)  {
        tmp_signals.pop_back();
        remapped_nodes[n][n] = -1;
      }
    } 
    if ( auto cand = associativity_candidate( ocs[0], ocs[1], ocs2[1], ocs2[2], ocs2[0] ); !mappable && cand ) {
      const auto& [x, y, z, u, assoc] = *cand;
      auto opt = ntk.create_maj( z, assoc ? u : x, ntk.create_maj( x, y, u ) );
     
      auto new_node = ntk.get_node(opt);
      tmp_signals.push_back(opt);
      remapped_nodes[n][n] = tmp_signals.size()-1; 

      if ( count_path_to_node (n, n, leaf, 0, true) == 1) {
        //std::cout << "\t\t\t\tassociativity 3 " << n << " with " << new_node << "\n";
        mappable = true;
        new_tt = flip_carry_user_tt(n, n, cuts.cuts(index)[0]);
        //ntk.set_cell_function(n, new_tt);
        //print_lut(n,n, false);
        //print_lut(n,n, true);
      }
      if (!mappable)  {
        tmp_signals.pop_back();
        remapped_nodes[n][n] = -1;
      }
    } 
 
    if(curr_tt != new_tt) {
      std::cout << "\t\t\t\t\tbefore: "; kitty::print_hex(curr_tt); std::cout <<"\n";
      std::cout << "\t\t\t\t\tafter: "; kitty::print_hex(new_tt); std::cout <<"\n";
      assert(0 && "After mapping, wrong truth table.\n");
    }


    return mappable;
  }

  using candidate_t = std::tuple<signal<Ntk>, signal<Ntk>, signal<Ntk>, signal<Ntk>, bool>;
  std::optional<candidate_t> associativity_candidate( signal<Ntk> const& v, signal<Ntk> const& w, signal<Ntk> const& x, signal<Ntk> const& y, signal<Ntk> const& z ) const
  {
    if ( v.index == x.index )
    {
      return candidate_t{w, y, z, v, v.complement == x.complement};
    }
    if ( v.index == y.index )
    {
      return candidate_t{w, x, z, v, v.complement == y.complement};
    }
    if ( w.index == x.index )
    {
      return candidate_t{v, y, z, w, w.complement == x.complement};
    }
    if ( w.index == y.index )
    {
      return candidate_t{v, x, z, w, w.complement == y.complement};
    }

    return std::nullopt;
  }


  bool distributivity( node<Ntk> const& n, uint32_t leaf, uint32_t rand ) {

    // Get order of children (depends on rand)
    const auto ocs = ordered_children( n, rand );

    // Check that the 3rd input has children to consider
    if ( !ntk.is_maj( ntk.get_node( ocs[2] ) ) || !ntk.is_maj( ntk.get_node( ocs[1] ) ) )
      return false;

    auto index = ntk.node_to_index(n);
    for ( auto cut_index: cuts.cuts(index)[0] ) {
      if (ntk.get_node(ocs[2]) == cut_index) return false;
      if (ntk.get_node(ocs[1]) == cut_index) return false;
    } 
    // Get children of 3rd input
    auto ocs1 = ordered_children( ntk.get_node( ocs[1] ), 0 );
    auto ocs2 = ordered_children( ntk.get_node( ocs[2] ), 0 );

    // Propagate inverter?
    if ( ntk.is_complemented( ocs[2] ) ) {
      ocs2[0] = !ocs2[0];
      ocs2[1] = !ocs2[1];
      ocs2[2] = !ocs2[2];
    }
    if ( ntk.is_complemented( ocs[1] ) ) {
      ocs1[0] = !ocs1[0];
      ocs1[1] = !ocs1[1];
      ocs1[2] = !ocs1[2];
    }

    auto curr_tt = flip_carry_user_tt(n,n, cuts.cuts(index)[0]);
    auto new_tt = curr_tt;
    bool mappable = false;
    
    if ( auto cand = distributivity_candidate(ocs1[0], ocs1[1], ocs1[2], ocs2[0], ocs2[1], ocs2[2]); cand ) {
      const auto& [x, y, u, v] = *cand;
      auto opt = ntk.create_maj( x, y, ntk.create_maj( u, v, ocs[0] ) );

      auto new_node = ntk.get_node(opt);
      tmp_signals.push_back(opt);
      remapped_nodes[n][n] = tmp_signals.size()-1; 

      //std::cout << "\t\t\t\tdistributivity " << n << " with " << new_node << "\n";
      if ( count_path_to_node (n, n, leaf, 0, true) == 1) {
        mappable = true;
        new_tt = flip_carry_user_tt(n, n, cuts.cuts(index)[0]);
        //ntk.set_cell_function(n, new_tt);
      }
      if (!mappable)  {
        tmp_signals.pop_back();
        remapped_nodes[n][n] = -1;
      }
    } 

    if(curr_tt != new_tt) {
      std::cout << "\t\t\t\t\tbefore: "; kitty::print_hex(curr_tt); std::cout <<"\n";
      std::cout << "\t\t\t\t\tafter: "; kitty::print_hex(new_tt); std::cout <<"\n";
      assert(0 && "After mapping, wrong truth table.\n");
    }
    return mappable;
  }

  using candidate2_t = std::tuple<signal<Ntk>, signal<Ntk>, signal<Ntk>, signal<Ntk>>;
  std::optional<candidate2_t> distributivity_candidate ( signal<Ntk> const& u, signal<Ntk> const& v, signal<Ntk> const& w, signal<Ntk> const& a, signal<Ntk> const& b, signal<Ntk> const& c) const
  {
    if ( a == u && b == v ) {
      return candidate2_t{a, b, c, w};
    }
    if ( a == u && b == w ) {
      return candidate2_t{a, b, c, v};
    }
    if ( a == v && b == w ) {
      return candidate2_t{a, b, c, u};
    }
    if ( a == v && b == u ) {
      return candidate2_t{a, b, c, w};
    }
    if ( a == w && b == u ) {
      return candidate2_t{a, b, c, v};
    }
    if ( a == w && b == v ) {
      return candidate2_t{a, b, c, u};
    }

    if ( a == u && c == v ) {
      return candidate2_t{a, c, b, w};
    }
    if ( a == u && c == w ) {
      return candidate2_t{a, c, b, v};
    }
    if ( a == v && c == w ) {
      return candidate2_t{a, c, b, u};
    }
    if ( a == v && c == u ) {
      return candidate2_t{a, c, b, w};
    }
    if ( a == w && c == u ) {
      return candidate2_t{a, c, b, v};
    }
    if ( a == w && c == v ) {
      return candidate2_t{a, c, b, u};
    }

    if ( b == u && c == v ) {
      return candidate2_t{b, c, a, w};
    }
    if ( b == u && c == w ) {
      return candidate2_t{b, c, a, v};
    }
    if ( b == v && c == w ) {
      return candidate2_t{b, c, a, u};
    }
    if ( b == v && c == u ) {
      return candidate2_t{b, c, a, w};
    }
    if ( b == w && c == u ) {
      return candidate2_t{b, c, a, v};
    }
    if ( b == w && c == v ) {
      return candidate2_t{b, c, a, u};
    }

    return std::nullopt;
  }


  std::array<signal<Ntk>, 3> ordered_children( node<Ntk> const& n, uint32_t rand ) const
  {
    std::array<signal<Ntk>, 3> children;

    ntk.foreach_fanin( n, [&children, rand]( auto const& f, auto i ) {
      if (rand == 0) {
        if (i==0) children[0] = f;
        else if (i == 1) children[1] = f;
        else children[2] = f;
      } else if (rand == 1) {
        if (i==0) children[0] = f;
        else if (i == 1) children[2] = f;
        else children[1] = f;
      } else if (rand ==2) {
        if (i==0) children[1] = f;
        else if (i == 1) children[2] = f;
        else children[0] = f;
      } else if (rand ==3) {
        if (i==0) children[1] = f;
        else if (i == 1) children[0] = f;
        else children[2] = f;
      } else if (rand == 4) {
        if (i==0) children[2] = f;
        else if (i == 1) children[0] = f;
        else children[1] = f;
      } else {
        if (i==0) children[2] = f;
        else if (i == 1) children[1] = f;
        else children[0] = f;
      } 
    } );

    //std::sort( children.begin(), children.end(), [this]( auto const& c1, auto const& c2 ) {
    //  return ntk.level( ntk.get_node( c1 ) ) < ntk.level( ntk.get_node( c2 ) );
    //} );
    return children;
  }

  bool try_other_graph ( uint32_t index, uint32_t leaf ) {
    uint32_t n = ntk.index_to_node(index);
    //std::cout << "\t\t\tRemapping " << index << " with " << remapped_nodes[n][n] << "\n";
    if (remapped_nodes[n][n] != -1) return true;

    if (relevance(n, leaf)) {
      //std::cout << "\t\t\t\t->Success\n";
      //print_lut(n,n, true);
      return true;
    }
    //std::cout << "\tdone relevance " << remapped_nodes[n][n] << "\n";
    for (uint32_t i = 0; i < 6; i++) {
      if (associativity(n, leaf, i)) {
        carry_nodes.resize(ntk.size());
        carry_cut_list.resize(ntk.size());
        delays.resize(ntk.size());
        //std::cout << "\t\t\t\t->Success\n";
        //print_lut(n,n, true);
        return true;
      } else if (distributivity (n, leaf, i)) {
        carry_nodes.resize(ntk.size());
        carry_cut_list.resize(ntk.size());
        delays.resize(ntk.size());
        //std::cout << "\t\t\t\t->Success\n";
        //print_lut(n,n, true);
        return true;
      }
    }
    //std::cout << "\tdone all" << remapped_nodes[n][n] << "\n";
    //std::cout << "\t\t\t\thow to update?\n";
    //print_lut(n,n,false);
    return false;
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
      for (uint32_t i = 0; i < carry_path.size(); i++) {
        std::cout << carry_path[i] << "(" << carry_driver_nodes[carry_path[i]] << "),";
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
        total_critical_path_count += critical_path_count;
      }
    });
    std::cout << "Total critical path is " << total_critical_path_count << " out of " <<  total_path_count << "\n";
    return total_critical_path_count;
  }

  void print_cut_cost (uint32_t index, uint32_t i1, uint32_t i2, uint32_t icut1, uint32_t icut2, double cost) {
    std::cout << "\tCut for index " << index << " "<< icut1 << " " << icut2 <<": ";
    std::cout << cost << "\n";
    std::cout << "\t\t" << i1 << ":" << cuts.cuts(i1)[icut1] << " " << cut_delay(cuts.cuts(i1)[icut1]) + LUT_ADDER_DELAY << "\n";
    std::cout << "\t\t" << i2 << ":" << cuts.cuts(i2)[icut2] << " " << cut_delay(cuts.cuts(i2)[icut2]) + LUT_ADDER_DELAY << "\n";
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
      worst_delay = delays[carry_driver_nodes[index]]+CARRY_DELAY;
      critical_index = carry_driver_nodes[index];
      for ( auto leaf : carry_cut_list[index] ) {
        std::cout << leaf;
        if (carry_driver_nodes[index] != leaf && delays[leaf]+LUT_ADDER_DELAY > worst_delay) {
          critical_index = leaf;
          worst_delay = delays[leaf]+LUT_ADDER_DELAY;
          std::cout << "*";
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
      //std::cout << "PO " << index << " has delay " << delays[index] << "\n";
      uint32_t curr_delay = delays[index];
      // THIS IS FOR THE ADDITIONAL INVERTER BUT NOT INCLUDING IT
      //if (is_a_carry_node(ntk.get_node(s))) curr_delay += LUT_DELAY;
      if (curr_delay > worst_delay) {
        critical_index = index;
        worst_delay = curr_delay;
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

  void print_params() {
    std::cout << "Mapping params:\n";
    std::cout << "\txilinx_arch " << ps.xilinx_arch << "\n";
    std::cout << "\tdelay " << ps.delay << "\n";
    std::cout << "\tcarry_mapping " << ps.carry_mapping << "\n";
    std::cout << "\tlength " << ps.length << "\n";
    std::cout << "\tdirect " << ps.direct << "\n";
    std::cout << "\tnum_crit_path " << ps.num_crit_path << "\n";

    
    /*for ( auto const& n : top_order )
    {
      if ( ntk.is_constant( n ) || ntk.is_pi( n ) || is_a_carry_node( n ) )
        continue;
      uint32_t cut_size = cuts.cuts(ntk.node_to_index(n)).size();
      if (cut_size == 100) std::cout << n << " > 100 \n";
    }*/
  }

  void print_state()
  {
    if (ps.verbose && ps.verbosity > 4 ) {
      for ( auto i = 0u; i < ntk.size(); ++i )
      {
        std::cout << fmt::format( "*** Obj = {:>3} (node = {:>3})  FlowRefs = {:5.2f}  MapRefs = {:>2}  Flow = {:5.2f}  Delay = {:>3} Carry = {}\n", \
          i, ntk.index_to_node( i ), flow_refs[i], map_refs[i], flows[i], delays[i], carry_nodes[i]);
      }
    }
    std::cout << fmt::format( "PS: Level = {}  Area = {}\n", delay/LUT_DELAY, area );
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
  std::vector<uint32_t> num_critical_path_through_node;
  std::vector<uint32_t> num_critical_map_refs;
  std::vector<bool> direct_input_nodes;
  std::vector<std::vector<int>> remapped_nodes;
  std::vector<bool> remapped_nodes_flip;
  std::vector<signal<Ntk>> tmp_signals;
  std::uint32_t tmp_size;

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


