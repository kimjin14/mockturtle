/* mockturtle: C++ logic network library
 * Copyright (C) 2018  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/*!
  \file carry_depth_view.hpp
  \brief Implements depth and level for a network

  \author Mathias Soeken
*/

#pragma once

#include <cstdint>
#include <vector>

#include "../traits.hpp"
#include "../utils/node_map.hpp"
#include "immutable_view.hpp"

namespace mockturtle
{

/*! \brief Implements `depth` and `level` methods for networks.
 *
 * This view computes the level of each node and also the depth of
 * the network.  It implements the network interface methods
 * `level` and `depth`.  The levels are computed at construction
 * and can be recomputed by calling the `update_levels` method.
 *
 * **Required network functions:**
 * - `size`
 * - `get_node`
 * - `visited`
 * - `set_visited`
 * - `foreach_fanin`
 * - `foreach_po`
 *
 * Example
 *
   \verbatim embed:rst

   .. code-block:: c++

      // create network somehow
      aig_network aig = ...;

      // create a depth view on the network
      carry_depth_view aig_depth{aig};

      // print depth
      std::cout << "Depth: " << aig_depth.depth() << "\n";
   \endverbatim
 */
template<typename Ntk, bool has_depth_interface = has_depth_v<Ntk>&& has_level_v<Ntk>&& has_update_levels_v<Ntk>>
class carry_depth_view
{
};

template<typename Ntk>
class carry_depth_view<Ntk, true> : public Ntk
{
public:
  carry_depth_view( Ntk const& ntk ) : Ntk( ntk )
  {
  }
};

template<typename Ntk>
class carry_depth_view<Ntk, false> : public Ntk
{
public:
  using storage = typename Ntk::storage;
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;

  /*! \brief Standard constructor.
   *
   * \param ntk Base network
   * \param count_complements Count inverters as 1
   */
  explicit carry_depth_view( Ntk const& ntk, bool count_complements = false )
      : Ntk( ntk ),
        _count_complements( count_complements ),
        _levels( ntk ),
        _levels_type( ntk )
  {
    static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
    static_assert( has_size_v<Ntk>, "Ntk does not implement the size method" );
    static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
    static_assert( has_is_complemented_v<Ntk>, "Ntk does not implement the is_complemented method" );
    static_assert( has_visited_v<Ntk>, "Ntk does not implement the visited method" );
    static_assert( has_set_visited_v<Ntk>, "Ntk does not implement the set_visited method" );
    static_assert( has_foreach_po_v<Ntk>, "Ntk does not implement the foreach_po method" );
    static_assert( has_foreach_fanin_v<Ntk>, "Ntk does not implement the foreach_fanin method" );

    update_levels();
  }

  uint32_t depth() const
  {
    return _depth;
  }

  uint32_t level( node const& n ) const
  {
    return _levels[n];
  }

  void print_levels() const
  {

    uint32_t n_out_accum = 0;
    uint32_t n_out = 0;
    for ( uint32_t d = 0; d < (_depth/20)+1; d++)
    {
      n_out = 0;
      this->foreach_po( [&]( auto const& f ) {
        if ( _levels[this->get_node(f)]/20 == d )
          n_out++;
      });
      n_out_accum += n_out;

      std::cout << "Level =\t" << d << ".\t";
      std::cout << "COs =\t" << n_out << ".\t";
      std::cout << "\t" << std::fixed << std::setprecision(0) << (float)n_out_accum/this->num_pos()*100 << "\%.\n";
    } 

    n_out = 0;
    this->foreach_po( [&]( auto const& f ) {
      if ( _levels[this->get_node(f)]/20 == (_depth/20) )
        n_out++;
    });
    std::cout << "Worst LUT level " << _depth/20 << ", COs = " << n_out << "/" << this->num_pos() << "\n";
  }

  void set_level( node const& n, uint32_t level )
  {
    _levels[n] = level;
  }

  void update_levels()
  {
    _levels.reset( 0 );

    this->incr_trav_id();
    compute_levels();
  }

  void resize_levels()
  {
    _levels.resize();
  }

private:
  uint32_t compute_levels( node const& n )
  {

    if ( this->visited( n ) == this->trav_id() )
    {
      return _levels[n];
    }
    this->set_visited( n, this->trav_id() );

    if ( this->is_constant( n ) || this->is_pi( n ) )
    {
      return _levels[n] = 0;
    }

    uint32_t level{0};
    node slowest_node = n;

    this->foreach_fanin( n, [&]( auto const& f ) {
    
      auto clevel = compute_levels( this->get_node( f ) );
     
      // If you want to count inversion as a delay 
      if ( _count_complements && this->is_complemented( f ) )
      {
        clevel+=LUT_DELAY;
      }
     
      if (this->is_carry(n) && (this->get_carry_driver(n) == this->get_node(f))) {
        clevel += CARRY_DELAY;
      } else if (this->is_carry(n)) {
        clevel += LUT_ADDER_DELAY;
      } else {
        clevel += LUT_DELAY;
      }

      if (clevel > level) {
        level = clevel;
        slowest_node = this->get_node(f);
      }
      level = std::max( level, clevel );
    } );

    // Return the correct delay
    // From your leaves, you have the worst delay node. 
    // You look at how much delay your node will have from the worst delay node.
    // If you are a carry and the node came from carry, you would add carry delay
    // If you are a carry and the node came from a LUT, you would add the large LUT delay
    // If you are a lut, you would add regular delay 
    return _levels[n] = level;
  }

  void compute_levels()
  {
    _depth = 0;
    this->foreach_po( [&]( auto const& f ) {
      auto clevel = compute_levels( this->get_node( f ) );
      if ( _count_complements && this->is_complemented( f ) )
      {
        clevel++;
      }
      _depth = std::max( _depth, clevel );
    } );
    print_depth_of_all_nodes();
    print_critical_path();
    
  }

  void print_critical_path_helper (uint32_t index, uint32_t& lut_count, uint32_t& lut_carry_count, uint32_t& carry_count) {
    
    auto n = this->index_to_node(index);

    //if (this->is_carry(n) ) {
    //  std::cout << " -> " << index << "*(" << _levels[index]/LUT_DELAY << ")\n";
    //} else {
    //  std::cout << " -> " << index << "(" << _levels[index]/LUT_DELAY << ")\n";
    //}

    if (this->is_pi(n) || this->is_constant(n))
      return;
   
    uint32_t worst_delay = 0;
    auto worst_node = 0;
    uint32_t critical_index = 0;

    this->foreach_fanin( n, [&]( auto const& f ) {
        const auto leaf = this->node_to_index (this->get_node(f));
        if (_levels[leaf] > worst_delay) {
          critical_index = leaf;
          worst_delay = _levels[leaf];
          worst_node = this->get_node(f);
        }
    } );

    if (this->is_carry(n) && (this->get_carry_driver(n) == worst_node)) {
      std::cout << " -> " << index << "*(" << _levels[index]/LUT_DELAY << ")\n";
      carry_count++;
    } else if (this->is_carry(n)) {
      std::cout << " -> " << index << "*(" << _levels[index]/LUT_DELAY << ")\n";
      lut_carry_count++;
    } else {
      std::cout << " -> " << index << "(" << _levels[index]/LUT_DELAY << ")\n";
      lut_count++;
    }

    print_critical_path_helper (critical_index, lut_count, lut_carry_count, carry_count);

  }

  void print_depth_of_all_nodes () {

    this->foreach_node ( [&]( auto n ) {
      if ( this->is_constant( n ) || this->is_pi( n ))
        return;

      else if (this->is_carry(n))
        std::cout << "Node* " << n << "(" << _levels[n] << "): ";
      else 
        std::cout << "Node " << n << "(" << _levels[n] << "): ";
      this->foreach_fanin( n, [&]( auto fanin ) {
          std::cout << fanin << "(" << _levels[fanin] << ") ";
      } );
      std::cout << "\n";
    } );
  }

  // Print the critical MIG path
  void print_critical_path () {

    std::cout << "Critical path:\n";
  
    uint32_t lut_count = 0;
    uint32_t carry_count = 0;  
    uint32_t lut_carry_count = 0;  

    uint32_t worst_delay = 0;
    uint32_t critical_index = 0;
    this->foreach_po( [&]( auto s ) {
      const auto index = this->node_to_index( this->get_node( s ) );
      if (_levels[index] > worst_delay) {
        critical_index = index;
        worst_delay = _levels[index];
      }
    });
    print_critical_path_helper (critical_index, lut_count, lut_carry_count, carry_count);
    std::cout << "total: " << lut_count << " " << lut_carry_count << " " << carry_count << "\n";
  }

  bool _count_complements{false};
  node_map<uint32_t, Ntk> _levels;
  node_map<uint32_t, Ntk> _levels_type;
  uint32_t _depth;
};

template<class T>
carry_depth_view( T const& ) -> carry_depth_view<T>;

template<class T>
carry_depth_view( T const&, bool ) -> carry_depth_view<T>;

} // namespace mockturtle
