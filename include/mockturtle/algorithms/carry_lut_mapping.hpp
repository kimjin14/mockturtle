/*
  carry_lut_mapping.hpp

  MIG nodes annotated on "value" (data[0].h2)
    if value is 1, mapped to carry
    if value is 2, mapped to LUT before carry     
*/

#pragma once

#include <cstdint>

#include <fmt/format.h>

#include "../utils/stopwatch.hpp"
#include "../views/topo_view.hpp"
#include "cut_enumeration.hpp"
#include "cut_enumeration/mf_cut.hpp"

namespace mockturtle
{

/*! \brief Parameters for lut_mapping.
 *
 * The data structure `lut_mapping_params` holds configurable parameters
 * with default arguments for `lut_mapping`.
 */
struct lut_mapping_params
{
  lut_mapping_params()
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
  bool verbose{false};
};

/*! \brief Statistics for lut_mapping.
 *
 * The data structure `lut_mapping_stats` provides data collected by running
 * `lut_mapping`.
 */
struct lut_mapping_stats
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
struct lut_mapping_update_cuts
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
class lut_mapping_impl
{
public:
  using network_cuts_t = network_cuts<Ntk, StoreFunction, CutData>;
  using cut_t = typename network_cuts_t::cut_t;

public:
  lut_mapping_impl( Ntk& ntk, lut_mapping_params const& ps, lut_mapping_stats& st )
      : ntk( ntk ),
        ps( ps ),
        st( st ),
        flow_refs( ntk.size() ),
        map_refs( ntk.size(), 0 ),
        flows( ntk.size() ),
        delays( ntk.size() ),
        carry_lut_nodes( ntk.size(), 0 ),
        cuts( cut_enumeration<Ntk, StoreFunction, CutData>( ntk, ps.cut_enumeration_ps ) )
  {
    lut_mapping_update_cuts<CutData>().apply( cuts, ntk );
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

/*
    int count_total_node = 0;
    int count_aig_node = 0;
    int count_mig_node = 0;

    ntk.foreach_node ([&](auto n) {
      if (ntk.is_pi(n) || ntk.is_constant(n)) return;

      count_total_node++;
      bool mig = true;
      ntk.foreach_fanin(n, [&](auto c) {
         if (ntk.is_constant(ntk.get_node(c))) mig = false;
      });

      if (mig) count_mig_node++;
      else count_aig_node++;

    });

    std::cout << count_total_node << "," << count_aig_node << ","<< count_mig_node << "\n";
*/
    //print_state();

    // Currently only finding one longest path (max of 30 nodes).
    // Using critical_path vector
    //std::cout << "Finding paths and mapping carry.\n";
    //find_critical_paths();
    //init_carry_nodes();

    // Map pre carry logic to LUTs
    // Using carry_lut_nodes vector 
    //std::cout << "Mapping to LUTs before carry nodes.\n";
    //compute_carry_mapping();

    // TODO: still needs to figure out if LUT nodes are used elsewhere
    // it currently only checks if its used in another carry node
    // This may work... still need to verify
    //for ( auto const& n : top_order ) {
    //  if (!is_a_carry_node(n)) {
    //    ntk.foreach_fanin(n, [&](auto const& c) {
    //      auto index = ntk.node_to_index(ntk.get_node(c));
    //      carry_lut_nodes[index]--; 
    //    });
    //  }
    //}
    /*for (int i = 0; i < carry_lut_nodes.size(); i++) {
      std::cout << carry_lut_nodes[i] << " ";
    }
    std::cout << "\n";*/

    //std::cout << "LUT mapping starts.\n";
    set_mapping_refs<false>();

    while ( iteration < ps.rounds )
    {
      compute_mapping<false>();
    }

    while ( iteration < ps.rounds + ps.rounds_ela )
    {
      compute_mapping<true>();
    }

    //std::cout << "Mapping derivation starts.\n";
    
    // TODO: mapping derivation needs to account for carry
    derive_mapping();
  }

private:
  uint32_t cut_area( cut_t const& cut ) const
  {
    return static_cast<uint32_t>( cut->data.cost );
  }

  bool get_path(node<Ntk> n, uint32_t depth, uint32_t curr_depth) {

    if ( ntk.is_constant( n ) || ntk.is_pi( n ) ) {
      if (curr_depth == depth || curr_depth == 40) {
        carry_nodes.push_back(n);
        return true;
      } else return false;
    }

    bool longest = false;

 
    for (uint32_t i = 0; i < ntk.fanin_size(n); i++) {
      auto nchild = ntk.get_children(n,i);
      longest = get_path(nchild, depth, curr_depth+1);
      if (longest) {
        carry_nodes.push_back(n);
        break;
      }
    }
    return longest;
  }

  void find_critical_paths() {

    auto depth = depth_view<Ntk>(ntk).depth();
    std::cout << "\tFinding path with depth of " << depth << ".\n";

    for (uint32_t i = 0; i < ntk.num_pos(); i++) {
      auto n = ntk.get_po(i);
      if (get_path(n,depth,0)) {
        break;
      }
    }
  }

  void init_nodes()
  {
    ntk.foreach_node( [this]( auto n, auto ) {
      const auto index = ntk.node_to_index( n );

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
  void init_carry_nodes() {
    for (auto n_carry: carry_nodes) {
      if (!ntk.is_constant(n_carry) && !ntk.is_pi(n_carry)) {
        const auto index = ntk.node_to_index( n_carry );
        delays[index] = 0;
      }
    }
  }

  // Decide which nodes can fit into LUT
  // Set mapping ref for those nodes
  void compute_carry_mapping()
  {
    uint32_t special_map = 0;

    // Map nodes following the critical path in pairs (due to ALM)
    for (uint32_t i = 1; i < carry_nodes.size()-1; i+=2) {
      auto const& n_carryin = carry_nodes[i-1];
      auto const& n_first = carry_nodes[i];
      auto const& n_second = carry_nodes[i+1];
       
      compute_best_cut_carry( ntk.node_to_index(n_first), ntk.node_to_index(n_second), \
        ntk.node_to_index(n_carryin), &special_map);
    }
    if (carry_nodes.size()%2==0) {
      auto const& n_carryin = carry_nodes[carry_nodes.size()-1-1];
      auto const& n_first = carry_nodes[carry_nodes.size()-1];
      compute_best_cut_carry( ntk.node_to_index(n_first), 0, \
        ntk.node_to_index(n_carryin), &special_map);
    }
    std::cout << "\tNumber of special mapping of nodes: " << special_map << "\n";
  }

  bool is_a_carry_node(node<Ntk> n) {
    for (uint32_t j = 0; j < carry_nodes.size(); j++) {
      if(n == carry_nodes[j]) {
        return true;
      }
    }
    return false; 
  }

  bool is_in_carry_lut(node<Ntk> n) {
    if (carry_lut_nodes[ntk.node_to_index(n)]>0) return true;
    return false;
    /*for (uint32_t j = 0; j < carry_lut_nodes.size(); j++) {
      if(== carry_lut_nodes[j]) {
        return true;
      }
    }
    return false; 
    */
  }

  template<bool ELA>
  void compute_mapping()
  {
    for ( auto const& n : top_order )
    {
      if ( ntk.is_constant( n ) || ntk.is_pi( n ) || is_a_carry_node(n))
        continue;
      compute_best_cut<ELA>( ntk.node_to_index( n ) );
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
      if (!is_a_carry_node(ntk.get_node(s))) {
        delay = std::max( delay, delays[index] );
        if constexpr ( !ELA )
        {
          map_refs[index]++;
        }
      }
    } );

    // increase the map_refs for those not mapped to carry or carry lut
    for (auto const n: carry_nodes) {
      if (ntk.is_pi(n)) continue;
      ntk.foreach_fanin(n, [&](auto const& c) {
        // map if not a PI, not in carry lut, and not a carry node
        if (!ntk.is_pi(ntk.get_node(c)) && !is_in_carry_lut(ntk.get_node(c)) && !is_a_carry_node(ntk.get_node(c))) {
          const auto index = ntk.node_to_index(ntk.get_node(c));
          //std::cout << "index of " << index << " set to map\n";
          map_refs[index]++;    
        }
      });
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
      //std::cout << "current " << index << "\n";
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

  bool insert_unique_input (node<Ntk> node, uint32_t* index_array, uint32_t* curr_index, uint32_t index1, uint32_t index2, 
      uint32_t carryin, uint32_t carryout, uint32_t max_num) {

    bool in_list = false;
    auto const& index = ntk.node_to_index(node);
    //if ((index != carryin) | (index != carryout) | (index != index1) | (index != index2)) {

      for (uint32_t i = 0; i < *curr_index; i++) {
        if (index_array[i] == index) 
          in_list = true;
      }
      
      std::cout << "Inserting " << index << " to " << *curr_index << " position.\n";
      // if there isn't 8 unique inputs yet, add 
      if (*curr_index < max_num && in_list == false) {
        index_array[*curr_index] = index;
        (*curr_index)++;
        return true;
      } else if (*curr_index >= max_num && in_list == false) {
        return false;
      } else return true; 

    //} else {

      // one of the special inputs, doesn't count towards unique input count
      //return true;
    //}
    return false;
  }
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
        //unique_array[*n_shared] = input_array[0][i+1];  
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
        //unique_array[index] = input_array[1][i+1];  
        index++;
      }
      match = false;
    } 

    return index;
  }

  // This function adds child's inputs to an array for input legality checking
  // 
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
          child_index[curr_LUT] = leaf_node;
          //std::cout << "input_array[" << curr_LUT << "][" << curr_input+1 << "] = " << leaf_node << ";\n";
          curr_input++;
        } else {
          // Check the children of the cut for inputs
          for (uint32_t i_fanin = 0; i_fanin < ntk.fanin_size(leaf_node); i_fanin++) { 
            node<Ntk> child_leaf_node = ntk.get_children (leaf_node, i_fanin);  
            input_array[curr_LUT][0]++;
            input_array[curr_LUT][curr_input+1] = child_leaf_node;
            //std::cout << "input_array[" << curr_LUT << "][" << curr_input+1 << "] = " << child_leaf_node << ";\n";
            curr_input++;
          }
          child_index[curr_LUT] = leaf_node;
        }
        curr_LUT++;
        curr_input = 0;
      }
    }
  }

  // Takes 2 indices and the carryin node value and determines whether its
  // children can fit into the LUTs preceding FA
  int check_child_node (uint32_t index1, uint32_t index2, uint32_t carryin) {

    // Array holding the nodes to the halfs of ALM as separate 4-LUT
    uint32_t index1_inputs[2][5]; 
    uint32_t index2_inputs[2][5]; 
    node<Ntk> index1_child[2] = {0};
    node<Ntk> index2_child[2] = {0};

    // Each node should have 3 children, where 1 is mapped to carry in
    // 2 inputs should use LUTs with above restrictions
    index1_inputs[0][0] = 0;
    index1_inputs[1][0] = 0;
    index2_inputs[0][0] = 0;
    index2_inputs[1][0] = 0;
    insert_child_to_array (index1_inputs, index1_child, index1, carryin);
    if (index2 != 0) {
      insert_child_to_array (index2_inputs, index2_child, index2, index1);
    }     
   /* 

    index1_inputs[0][0] = 4;
    index1_inputs[0][1] = 0; 
    index1_inputs[0][2] = 1; 
    index1_inputs[0][3] = 2; 
    index1_inputs[0][4] = 3; 

    index1_inputs[1][0] = 4;
    index1_inputs[1][1] = 0; 
    index1_inputs[1][2] = 1; 
    index1_inputs[1][3] = 2; 
    index1_inputs[1][4] = 4; 

    index2_inputs[0][0] = 4;
    index2_inputs[0][1] = 0; 
    index2_inputs[0][2] = 1; 
    index2_inputs[0][3] = 4; 
    index2_inputs[0][4] = 5; 

    index2_inputs[1][0] = 4;
    index2_inputs[1][1] = 0; 
    index2_inputs[1][2] = 1; 
    index2_inputs[1][3] = 4; 
    index2_inputs[1][4] = 3; 
*/
/*
    std::cout << "LUT0_0 ";
    for (uint32_t i = 0; i < index1_inputs[0][0]; i++) {
      //index1_inputs[0][i+1] = i;
      std::cout << index1_inputs[0][i+1] << " ";
    }
    std::cout << "\n";
    std::cout << "LUT0_1 ";
    for (uint32_t i = 0; i < index1_inputs[1][0]; i++) {
      //index1_inputs[1][i+1] = i;
      std::cout << index1_inputs[1][i+1] << " ";
    }
    std::cout << "\n";
    std::cout << "LUT1_0 ";
    for (uint32_t i = 0; i < index2_inputs[0][0]; i++) {
      //index2_inputs[0][i+1] = i+3;
      std::cout << index2_inputs[0][i+1] << " ";
    }
    std::cout << "\n";
    std::cout << "LUT1_1 ";
    for (uint32_t i = 0; i < index2_inputs[1][0]; i++) {
      //index2_inputs[1][i+1] = i+3;
      std::cout << index2_inputs[1][i+1] << " ";
    }
    std::cout << "\n";
*/
    // Check if inputs are legal
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

    
    // node 1 and node 2 both exceed unique input counts
    // nothing can map in this case
    if (n_unique_index1 > 5 && n_unique_index2 > 5) {
      return 0;
    }

    // at least one of the index has less than 5 unique inputs
    // if one of them has more it cannot map, the other can definitely map
    if (n_unique_index1 > 5) {
      //std::cout << "index1 node cannot fit\n";
      for (uint32_t i = 0; i < 2; i++) { 
        if (index2_child[i] != 0) {
          carry_lut_nodes[index2_child[i]]++;
          //carry_lut_nodes.push_back(index2_child[i]); 
        }
      }
      for (uint32_t i = 0; i < 2; i++) { 
        if (index1_child[i] != 0) {
          carry_lut_nodes[index1_child[i]]--;
          //carry_lut_nodes.push_back(index2_child[i]); 
        }
      }
      return 1;
    }

    // index 2 definitely cannot map but index 1 is able to
    if (index2 == 0 || n_unique_index2 > 5) {
      //std::cout << "index2 node cannot fit\n";
      for (uint32_t i = 0; i < 2; i++) { 
        if (index1_child[i] != 0) {
          carry_lut_nodes[index1_child[i]]++;
          //carry_lut_nodes.push_back(index1_child[i]); 
        }
      }
      for (uint32_t i = 0; i < 2; i++) { 
        if (index2_child[i] != 0) {
          carry_lut_nodes[index2_child[i]]--;
          //carry_lut_nodes.push_back(index2_child[i]); 
        }
      }
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
/*
    for (int i = 0; i < 5; i++) {
      std::cout << abce0f0[i] << " ";
    }
    std::cout << "\n";
    for (int i = 0; i < 5; i++) {
      std::cout << abde1f1[i] << " ";
    }
    std::cout << "\n";
    std::cout << n_unique_index1 << "\n";
    std::cout << n_unique_index2 << "\n";
    std::cout << n_matches << "\n";
*/
    if (n_unique_index1 + n_unique_index2 - n_matches <= 8) {
    //if (((n_unique_index1 == 5 || n_unique_index2 == 5) && n_matches >= 2) || 
    //    ((n_unique_index1 == 4 || n_unique_index2 == 4) && n_matches >= 1) || 
    //    (n_unique_index1 <=3 && n_unique_index2 <=3) ) {
      for (uint32_t i = 0; i < 2; i++) { 
        if (index1_child[i] != 0) {
          carry_lut_nodes[index1_child[i]];
          //carry_lut_nodes.push_back(index1_child[i]); 
        }
        if (index2_child[i] != 0) {
          carry_lut_nodes[index2_child[i]];
          //carry_lut_nodes.push_back(index2_child[i]); 
        }
      }
      return 2;
    }
    std::cout << "Nothing fits\n";

    for (uint32_t i = 0; i < 2; i++) { 
      if (index1_child[i] != 0) {
        carry_lut_nodes[index1_child[i]]--;
        //carry_lut_nodes.push_back(index2_child[i]); 
      }
    }
    for (uint32_t i = 0; i < 2; i++) { 
      if (index2_child[i] != 0) {
        carry_lut_nodes[index2_child[i]]--;
        //carry_lut_nodes.push_back(index2_child[i]); 
      }
    }

    return 0;
  } 

  bool cut_check_legality( cut_t const& cut1, cut_t const& cut2, uint32_t index1, uint32_t index2, 
      uint32_t carryin, uint32_t carryout ) {

    // Checking the usage of (Stratix V) LUTs before mapping MIG node to carry
    // 1. inputs to the cuts (unless it is a PI or carry signal) cannot be more than 8
    // 2. each node can have at most 2 pairs of unique inputs

    // Each node can have 4 unique inputs
    uint32_t curr_inputs = 0;
    uint32_t unique_inputs[8];

    for (uint32_t leaf : cut1 ) {

      auto const& leaf_node = ntk.index_to_node(leaf);

      // If its one of carry connection, skip
      if (ntk.is_constant(leaf_node) | (leaf_node == carryin) | (leaf_node == carryout) | \
          (leaf_node == index1) | (leaf_node == index2)) 
        continue;

      // If leaf is input, no need to check children
      if (ntk.is_pi(leaf_node)) {
        if (!insert_unique_input (leaf_node, unique_inputs, &curr_inputs,
            index1, index2, carryin, carryout, 6)) return false;
      } else {

        // Check the children of the cut for inputs
        for (uint32_t i_fanin = 0; i_fanin < ntk.fanin_size(leaf_node); i_fanin++) { 
          node<Ntk> child_leaf_node = ntk.get_children (leaf_node, i_fanin);  
  
          if (!insert_unique_input (child_leaf_node, unique_inputs, &curr_inputs,
              index1, index2, carryin, carryout, 6)) return false;
        }
      }
    }
    std::cout << "\tPassed 4 unique input test 1\n";
    
    curr_inputs = 6;
    for (auto leaf : cut2 ) {

      // check the children of the cut
      auto const& leaf_node = ntk.index_to_node(leaf);
       
      // If in carry connection, skip
      if (ntk.is_constant(leaf_node) | (leaf_node == carryin) | (leaf_node == carryout) | \
          (leaf_node == index1) | (leaf_node == index2)) 
        continue;

      if (ntk.is_pi(leaf_node)) {
        if (!insert_unique_input (leaf_node, unique_inputs, &curr_inputs,
            index1, index2, carryin, carryout, 8)) return false;
      } else {

        // if leaf is input, no need to check children
        // check the children of the cut
        for (uint32_t i_fanin = 0; i_fanin < ntk.fanin_size(leaf_node); i_fanin++) { 
          //std::cout << "\t\t";
          node<Ntk> child_leaf_node = ntk.get_children (leaf_node, i_fanin);  
          //std::cout << child_leaf_node << "\n";
  
          if (!insert_unique_input (child_leaf_node, unique_inputs, &curr_inputs,
              index1, index2, carryin, carryout, 8)) return false;
        }
      }
    }
    std::cout << "\tPassed 4 unique input test 2\n";
    //std::cout << "\tPassed 8 shared input test\n";

    //for (uint32_t i = 0; i < curr_inputs; i++)
    //  std::cout << unique_inputs[i] << " " ;
    //std::cout << "\n";

    return true;
  }

  std::pair<float, uint32_t> cut_flow_carry( cut_t const& cut1, cut_t const& cut2)
  {
    uint32_t time{0u};
    float flow{0.0f};
    for ( auto leaf : cut1 )
    {
      time = std::max( time, delays[leaf] );
      flow += flows[leaf];
    }
    for ( auto leaf : cut2 )
    {
      time = std::max( time, delays[leaf] );
      flow += flows[leaf];
    }
    return {flow + cut_area( cut1 ) + cut_area( cut2), time + 1u};
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

    return {flow + cut_area( cut ), time + 1u};
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

  /* this computes the best cut for two nodes at the same time
    since both nodes will be packed to the same ALM */
  void compute_best_cut_carry( uint32_t index1, uint32_t index2, \
      uint32_t index_carryin, uint32_t* n_special_map)
  {
    uint32_t n_mapped = check_child_node(index1, index2, index_carryin);
    *n_special_map += n_mapped;
    return;
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
        //std::cout << "best_cut = " << best_cut << "\n";
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

  // Add to carry mapping, which tells the mapping view that this node is in carry
  // Add to regular mapping with inputs correctly representing either
  //    mapped child node inputs, itself (just a wire) 
  void map_paths_to_carry_chain()
  {

    for (auto const n: carry_nodes){
      if ( ntk.is_constant( n ) || ntk.is_pi( n )) 
        continue;
      //std::cout << "\tCarry Mapping " << n << "\n";

      
      std::vector<node<Ntk>> nodes;

      //std::cout << "\t\tAdding: ";
      for (uint32_t c = 0; c < ntk.fanin_size(n); c++)
      {
        auto nchild = ntk.get_children(n,c);
        //std::cout << nchild << "-> ";
        if (ntk.is_pi(nchild) ){
          nodes.push_back(nchild);
          //std::cout << nchild << " ";
        } else if (is_a_carry_node(nchild)) {
          continue;
        } else if (is_in_carry_lut(nchild)) { 
          //std::cout << "in carry lut " ;
          for (uint32_t i = 0; i < ntk.fanin_size(nchild); i++) {
            auto nchildinput = ntk.get_children(nchild,i);
            nodes.push_back(nchildinput); 
            //std::cout << nchildinput << " ";
          }
        } else {
          nodes.push_back(nchild);
          //std::cout << nchild << " ";
        }
      }

     //std::cout << "\n"; 

      //ntk.add_to_mapping (n, nodes.begin(), nodes.end()); 
      // 1 is for which carry chain this belongs to. For now, just 1
      ntk.add_to_mapping( n, nodes.begin(), nodes.end() );
      ntk.add_to_carry_mapping (n, 1); 
    }
    /*for (auto const n: carry_lut_nodes)
    {
      if (!is_in_carry_lut(n)) continue;
      
      std::vector<node<Ntk>> nodes;

      for (uint32_t c = 0; c < ntk.fanin_size(n); c++)
      {
        auto nchild = ntk.get_children(n,c);
        nodes.push_back(nchild);
      }

      ntk.add_to_mapping (n, nodes.begin(), nodes.end()); 
    }*/
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
      //if (printtt)
      //std::cout << n << " -> " << c << "\n";
      for ( auto const& l : bestcut ) {
        if (l == ntk.node_to_index(c)) {
          if (*max_depth < depth)
            *max_depth = depth;
          found = true;
        } 
      } 
      if (!found && !ntk.is_pi(c) && !ntk.is_constant(c) && !(ntk.visited(c)>0)) {
        //if (printtt)
        //std::cout << "\tn " << n << " c " << c << "\n";
        //if (printtt)
        //std::cout << c << ",";
        (*num_nodes)++;
        find_num_nodes(c, bestcut, num_nodes, depth+1, max_depth, printtt);
      }
    }
  }

  // n - node you are currently investigating
  // index - index of one of the cut I'm looking at 
  /*void find_depth(node<Ntk> n, uint32_t index, uint32_t depth, \
      uint32_t* max_depth, uint32_t* num_nodes) {
    for (uint32_t i = 0; i < ntk.fanin_size(n); i++) {
      auto c = ntk.get_children(n,i);
      if (index==ntk.node_to_index(c)) {
        if (*max_depth < depth) 
          *max_depth = depth;
        //std::cout << "\t" << index << "\n";
      } else {
        //(*num_nodes)++;
        //std::cout << "\t" << index << " "<< ntk.node_to_index(c) << "\n";
        find_depth(c, index, depth+1, max_depth, num_nodes);
      }
    }
  }*/


  void derive_mapping()
  {
    ntk.clear_mapping();

    map_paths_to_carry_chain();

    for ( auto const& n : top_order )
    {
      if ( ntk.is_constant( n ) || ntk.is_pi( n ) || ntk.is_cell_root( n ) ) 
        continue;

      const auto index = ntk.node_to_index( n );

      if (is_a_carry_node(n))
        assert(map_refs[index] == 0);

      if ( map_refs[index] == 0 )
        continue;

      std::cout << "\tRegular Mapping " << n << "\n";
      std::vector<node<Ntk>> nodes;
      for ( auto const& l : cuts.cuts( index ).best() ) {
        nodes.push_back( ntk.index_to_node( l ) );
      }
      uint32_t depth = 0;
      uint32_t num_nodes = 1;
      find_num_nodes (n, cuts.cuts( index ).best(), &num_nodes, 1, &depth, false);
      ntk.clear_visited();
      
      ntk.add_to_mapping( n, nodes.begin(), nodes.end() );

      if constexpr ( StoreFunction )
      {
        ntk.set_cell_function( n, cuts.truth_table( cuts.cuts( index ).best() ) );
      }
    }
  }

  void print_state()
  {
    for ( auto i = 0u; i < ntk.size(); ++i )
    {
      std::cout << fmt::format( "*** Obj = {:>3} (node = {:>3})  FlowRefs = {:5.2f}  MapRefs = {:>2}  Flow = {:5.2f}  Delay = {:>3}\n", i, ntk.index_to_node( i ), flow_refs[i], map_refs[i], flows[i], delays[i] );
    }
    std::cout << fmt::format( "Level = {}  Area = {}\n", delay, area );
  }

private:
  Ntk& ntk;
  lut_mapping_params const& ps;
  lut_mapping_stats& st;

  uint32_t iteration{0}; /* current mapping iteration */
  uint32_t delay{0};     /* current delay of the mapping */
  uint32_t area{0};      /* current area of the mapping */
  //bool ela{false};       /* compute exact area */

  std::vector<node<Ntk>> top_order;
  std::vector<float> flow_refs;
  std::vector<uint32_t> map_refs;
  std::vector<float> flows;
  std::vector<uint32_t> delays;
  std::vector<int> carry_lut_nodes;
  network_cuts_t cuts;

  std::vector<uint32_t> tmp_area; /* temporary vector to compute exact area */

  // Contains the list of nodes to be mapped to carry 
  std::vector<node<Ntk>> carry_nodes;
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
void carry_lut_mapping( Ntk& ntk, lut_mapping_params const& ps = {}, lut_mapping_stats* pst = nullptr )
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

  lut_mapping_stats st;
  detail::lut_mapping_impl<Ntk, StoreFunction, CutData> p( ntk, ps, st );
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