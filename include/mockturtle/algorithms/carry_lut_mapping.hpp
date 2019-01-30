/*
  carry_lut_mapping.hpp
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

    find_critical_paths();
    set_mapping_refs<false>();
    print_state();

    while ( iteration < ps.rounds )
    {
      compute_mapping<false>();
      //print_state();
    }

    while ( iteration < ps.rounds + ps.rounds_ela )
    {
      compute_mapping<true>();
      //print_state();
    }

    derive_mapping();
    
  }

private:
  uint32_t cut_area( cut_t const& cut ) const
  {
    return static_cast<uint32_t>( cut->data.cost );
  }

  bool get_path(node<Ntk> n, uint32_t depth, uint32_t curr_depth) {

    if ( ntk.is_constant( n ) || ntk.is_pi( n ) ) {
      if (curr_depth == depth) {
        //critical_path.push_back(n);
        return true;
      } else return false;
    }

    bool longest = false;
 
    for (uint32_t i = 0; i < ntk.fanin_size(n); i++) {
      auto nchild = ntk.get_children(n,i);
      longest = get_path(nchild, depth, curr_depth+1);
      if (longest) {
        critical_path.push_back(n);
        break;
      }
    }
    return longest;
  }

  // find longest path and place it on carry
  // currently: finds the first longest path
  void find_critical_paths() {

    auto depth = depth_view<Ntk>(ntk).depth();
    std::cout << "finding depth " << depth << "\n";

    for (uint32_t i = 0; i < ntk.num_pos(); i++) {
      auto n = ntk.get_po(i);
      //std::cout << n << ":"; 
      if (get_path(n,depth,0)) {
        break;
      }
      //std::cout << "\n";
    }

    for (uint i = 0; i < critical_path.size(); i++) 
      std::cout << critical_path[i] << " "; 
    std::cout << "\n";
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

  template<bool ELA>
  void compute_mapping()
  {
    bool carry_index = false;
    for ( auto const& n : top_order )
    {
      if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
        continue;

      // If node is in carry chain, remove illegal cuts
      for (uint32_t i = 1; i < critical_path.size(); i+=2) {
        if (n == critical_path[i]) {

          carry_index = true;
          auto const& n_connected = critical_path[i+1];
          auto const& index_carryin = critical_path[i-1];
          uint32_t index_carryout = 0;
          if (i+2 < critical_path.size())
            index_carryout = ntk.node_to_index (critical_path[i+2]);
          std::cout << n << "&" << n_connected << "\n";
          compute_best_cut_carry<ELA>( ntk.node_to_index( n ), ntk.node_to_index (n_connected), \
            ntk.node_to_index(index_carryin), index_carryout );
        }
      }
      if (!carry_index) {
        std::cout << n << "\n";
        compute_best_cut<ELA>( ntk.node_to_index( n ) );
      }
    }
    set_mapping_refs<ELA>();
    //print_state();
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

      if constexpr ( !ELA )
      {
        map_refs[index]++;
      }
    } );

    // Anything going into the carry chain should "seem" like a po
    for (auto const n : critical_path) {
      if (ntk.is_pi(n)) continue;
      const auto index = ntk.node_to_index (n); 
      if constexpr ( !ELA )
      {
        map_refs[index]++;
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
      uint32_t carryin, uint32_t carryout ) {

    bool in_list = false;
    auto const& index = ntk.node_to_index(node);
    if ((index != carryin) | (index != carryout) | (index != index1) | (index != index2)) {

      for (uint32_t i = 0; i < *curr_index; i++) {
        if (index_array[i] == index) 
          in_list = true;
      }
      
      // if there isn't 8 unique inputs yet, add 
      if (*curr_index < 8 && in_list == false) {
        index_array[*curr_index] = index;
        (*curr_index)++;
        return true;
      } else if (*curr_index >= 8 && in_list == false) {
        return false;
      } else return true; 

    } else {

      // one of the special inputs, doesn't count towards unique input count
      return true;
    }
    return false;
  }

  bool cut_check_legality( cut_t const& cut1, cut_t const& cut2, uint32_t index1, uint32_t index2, 
      uint32_t carryin, uint32_t carryout ) {

    std::cout << "Cut check:\n";

    // CHECK #1
    // Each cut should be size for 4 LUT
    //if (cut1.size() > 4 ) return false;
    //if (cut2.size() > 4 ) return false;
    //std::cout << "\tPassed 4 LUT test\n";
    std::cout << "\t" << cut1 << "\n";
    std::cout << "\t" << cut2 << "\n";
  
    // CHECK #2
    // Inputs to these nodes in the cut must not have more than 8 unique inputs
    uint32_t curr_inputs = 0;
    uint32_t unique_inputs[8];
    std::cout << "\tChecking cut 1\n";
    for (uint32_t leaf : cut1 ) {

      auto const& leaf_node = ntk.index_to_node(leaf);

      // If in carry connection, skip
      if ((leaf_node == carryin) | (leaf_node == carryout) | \
          (leaf_node == index1) | (leaf_node == index2)) 
        continue;

      if (ntk.is_pi(leaf_node)) {
        if (!insert_unique_input (leaf_node, unique_inputs, &curr_inputs,
            index1, index2, carryin, carryout)) return false;
      } else {

        // if leaf is input, no need to check children
        // check the children of the cut
        for (uint32_t i_fanin = 0; i_fanin < ntk.fanin_size(leaf_node); i_fanin++) { 
          //std::cout << "\t\t";
          node<Ntk> child_leaf_node = ntk.get_children (leaf_node, i_fanin);  
          //std::cout << child_leaf_node << "\n";
  
          if (!insert_unique_input (child_leaf_node, unique_inputs, &curr_inputs,
              index1, index2, carryin, carryout)) return false;
        }
      }
    }
    std::cout << "\tChecking cut 2\n";
    for (auto leaf : cut2 ) {

      // check the children of the cut
      auto const& leaf_node = ntk.index_to_node(leaf);
       
      // If in carry connection, skip
      if ((leaf_node == carryin) | (leaf_node == carryout) | \
          (leaf_node == index1) | (leaf_node == index2)) 
        continue;

      if (ntk.is_pi(leaf_node)) {
        if (!insert_unique_input (leaf_node, unique_inputs, &curr_inputs,
            index1, index2, carryin, carryout)) return false;
      } else {

        // if leaf is input, no need to check children
        // check the children of the cut
        for (uint32_t i_fanin = 0; i_fanin < ntk.fanin_size(leaf_node); i_fanin++) { 
          //std::cout << "\t\t";
          node<Ntk> child_leaf_node = ntk.get_children (leaf_node, i_fanin);  
          //std::cout << child_leaf_node << "\n";
  
          if (!insert_unique_input (child_leaf_node, unique_inputs, &curr_inputs,
              index1, index2, carryin, carryout)) return false;
        }
      }
    }
    std::cout << "\tPassed 8 shared input test\n";

    for (uint32_t i = 0; i < curr_inputs; i++)
      std::cout << unique_inputs[i] << " " ;
    std::cout << "\n";

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
  template<bool ELA>
  void compute_best_cut_carry( uint32_t index1, uint32_t index2, uint32_t index_carryin, uint32_t index_carryout )
  {
    constexpr auto mf_eps{0.005f};

    float flow;
    uint32_t time{0};
    int32_t best_cut1{-1};
    int32_t best_cut2{-1};
    float best_flow{std::numeric_limits<float>::max()};
    uint32_t best_time{std::numeric_limits<uint32_t>::max()};
    int32_t cut_index1{-1};
    int32_t cut_index2{-1};

    //if constexpr ( ELA )
    //{
    //  if ( map_refs[index] > 0 )
    //  {
    //    cut_deref( cuts.cuts( index1 )[0] );
    //  }
    //}

    for (auto* cut1: cuts.cuts(index1)) {
      ++cut_index1;
      if ( cut1->size() == 1 )
        continue;
      cut_index2 = 0;
      for (auto* cut2: cuts.cuts(index2)) {
        if ( cut2->size() == 1 )
          continue;
        ++cut_index2;
        if constexpr ( ELA )
        {
          flow = static_cast<float>( cut_area_estimation( *cut1 ) );
        }
        else
        {
          if (cut_check_legality(*cut1, *cut2 , index1, index2, index_carryin, index_carryout)) {
            std::cout << "Valid mapping for LUTs before carry chain.\n";
  
            std::tie( flow, time ) = cut_flow_carry( *cut1, *cut2);
          }
        }

        if ( best_cut1 == -1 || best_flow > flow + mf_eps || ( best_flow > flow - mf_eps && best_time > time ) )
        {
          best_cut1 = cut_index1;
          best_cut2 = cut_index2;
          best_flow = flow;
          best_time = time;
        }
      }
    }

    //if constexpr ( ELA )
    //{
    //  if ( map_refs[index] > 0 )
    //  {
    //    cut_ref( cuts.cuts( index )[best_cut] );
    //  }
    //}
    //else
    //{
      map_refs[index1] = 0;
      map_refs[index2] = 0;
    //}
    //if constexpr ( ELA )
    //{
    //  best_time = cut_flow( cuts.cuts( index )[best_cut] ).second;
    //}
    delays[index1] = best_time;
    flows[index1] = best_flow / flow_refs[index1];
    delays[index2] = best_time;
    flows[index2] = best_flow / flow_refs[index2];

    if ( best_cut1 != 0 )
    {
      cuts.cuts( index1 ).update_best( best_cut1 );
    }
    if ( best_cut2 != 0 )
    {
      cuts.cuts( index2 ).update_best( best_cut2 );
    }
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

  void map_paths_to_carry_chain()
  {

    for (auto const n: critical_path)
    {
      std::vector<node<Ntk>> nodes;

      for (uint32_t c = 0; c < ntk.fanin_size(n); c++)
      {
        auto nchild = ntk.get_children(n,c);
        nodes.push_back(nchild);
      }

      ntk.add_to_mapping (n, nodes.begin(), nodes.end()); 
    }
  }

  void derive_mapping()
  {
    ntk.clear_mapping();
    //map_paths_to_carry_chain();

    for ( auto const& n : top_order )
    {
      if ( ntk.is_constant( n ) || ntk.is_pi( n ) || ntk.is_cell_root( n ) ) 
        continue;

      const auto index = ntk.node_to_index( n );
      if ( map_refs[index] == 0 )
        continue;

      std::cout << "here " << n << "\n";
      std::vector<node<Ntk>> nodes;
      for ( auto const& l : cuts.cuts( index ).best() )
      {
        nodes.push_back( ntk.index_to_node( l ) );
      }
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
  network_cuts_t cuts;

  std::vector<uint32_t> tmp_area; /* temporary vector to compute exact area */

  // Contains the list of nodes in the critical_path
  std::vector<node<Ntk>> critical_path;
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
