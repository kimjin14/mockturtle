/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2019  EPFL
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
  \file mapping_view.hpp
  \brief Implements mapping methods to create mapped networks

  \author Mathias Soeken
*/

#pragma once

#include <memory>
#include <vector>
#include <iomanip>

#include "../networks/detail/foreach.hpp"
#include "../traits.hpp"
#include "../utils/truth_table_cache.hpp"
#include "../views/immutable_view.hpp"

#include <kitty/dynamic_truth_table.hpp>

namespace mockturtle
{

namespace detail
{

template<bool StoreFunction>
struct mapping_view_storage;

template<>
struct mapping_view_storage<true>
{
  std::vector<uint32_t> mappings;
  std::vector<uint32_t> carry_mappings;
  std::vector<uint32_t> carry_luta_mappings;
  std::vector<uint32_t> carry_lutb_mappings;
  std::vector<bool> carry_lut_mapped;
  std::vector<bool> carry_lut_mapping_complemented;
  uint32_t mapping_size{0};
  uint32_t carry_mapping_size{0};
  uint32_t carry_lut_mapping_size{0};
  std::vector<uint32_t> functions;
  truth_table_cache<kitty::dynamic_truth_table> cache;
};

template<>
struct mapping_view_storage<false>
{
  std::vector<uint32_t> mappings;
  // 0 if not in carry chain, number indicates part of which carry chain
  std::vector<uint32_t> carry_mappings;
  // which luts are actually absorbed into 
  std::vector<uint32_t> carry_luta_mappings;
  std::vector<uint32_t> carry_lutb_mappings;
  std::vector<bool> carry_lut_mapped;
  std::vector<bool> carry_lut_mapping_complemented;
  uint32_t mapping_size{0};
  uint32_t carry_mapping_size{0};
  uint32_t carry_lut_mapping_size{0};
};

} // namespace detail

template<typename Ntk, bool StoreFunction>
inline constexpr bool implements_mapping_interface_v = has_has_mapping_v<Ntk> && (!StoreFunction || has_cell_function_v<Ntk>);

/*! \brief Adds mapping API methods to network.
 *
 * This view adds methods of the mapping API methods to a network.  It always
 * adds the functions `has_mapping`, `is_cell_root`, `clear_mapping`,
 * `num_cells`, `add_to_mapping`, `remove_from_mapping`, and
 * `foreach_cell_fanin`.  If the template argument `StoreFunction` is set to
 * `true`, it also adds functions for `cell_function` and `set_cell_function`.
 * For the latter case, this view requires more memory to also store the cells'
 * truth tables.
 *
 * **Required network functions:**
 * - `size`
 * - `node_to_index`
 *
 * Example
 *
   \verbatim embed:rst

   .. code-block:: c++

      // create network somehow
      aig_network aig = ...;

      // in order to apply mapping, wrap network in mapping view
      mapping_view mapped_aig{aig};

      // call LUT mapping algorithm
      lut_mapping( mapped_aig );

      // nodes of aig and mapped_aig are the same
      aig.foreach_node( [&]( auto n ) {
        std::cout << n << " has mapping? " << mapped_aig.is_cell_root( n ) << "\n";
      } );
   \endverbatim
 */
template<typename Ntk, bool StoreFunction = false, bool has_mapping_interface = implements_mapping_interface_v<Ntk, StoreFunction>>
class mapping_view
{
};

template<typename Ntk, bool StoreFunction>
class mapping_view<Ntk, StoreFunction, true> : public Ntk
{
public:
  mapping_view( Ntk const& ntk ) : Ntk( ntk )
  {
  }
};

template<typename Ntk, bool StoreFunction>
class mapping_view<Ntk, StoreFunction, false> : public immutable_view<Ntk>
{
public:
  using storage = typename Ntk::storage;
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;

  /*! \brief Default constructor.
   *
   * Constructs mapping view on another network.
   */
  mapping_view( Ntk const& ntk )
      : immutable_view<Ntk>( ntk ),
        _mapping_storage( std::make_shared<typename decltype( _mapping_storage )::element_type>() )
  {
    static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
    static_assert( has_size_v<Ntk>, "Ntk does not implement the size method" );
    static_assert( has_node_to_index_v<Ntk>, "Ntk does not implement the node_to_index method" );

    _mapping_storage->mappings.resize( ntk.size(), 0 );
    _mapping_storage->carry_mappings.resize( ntk.size(), 0 );
    _mapping_storage->carry_luta_mappings.resize( ntk.size(), 0 );
    _mapping_storage->carry_lutb_mappings.resize( ntk.size(), 0 );
    _mapping_storage->carry_lut_mapped.resize( ntk.size(), false );
    _mapping_storage->carry_lut_mapping_complemented.resize( ntk.size(), false );

    if constexpr ( StoreFunction )
    {
      /* insert 0 truth table */
      _mapping_storage->cache.insert( kitty::dynamic_truth_table( 0 ) );

      /* default each truth table to 0 */
      _mapping_storage->functions.resize( ntk.size(), 0 );
    }
  }

  void resize_mapping (uint32_t add) {

    if constexpr ( StoreFunction )
    {
      /* insert 0 truth table */
      _mapping_storage->cache.insert( kitty::dynamic_truth_table( 0 ) );

      /* default each truth table to 0 */
      _mapping_storage->functions.resize( this->size() + add, 0 );
    }
    _mapping_storage->mappings.resize( this->size() + add, 0);
    _mapping_storage->carry_mappings.resize( this->size() + add, 0);
    _mapping_storage->carry_luta_mappings.resize( this->size() + add, 0 );
    _mapping_storage->carry_lutb_mappings.resize( this->size() + add, 0 );
    _mapping_storage->carry_lut_mapped.resize( this->size() + add, false );
    _mapping_storage->carry_lut_mapping_complemented.resize( this->size() + add, false );
  }
  
  bool has_mapping() const
  {
    return _mapping_storage->mapping_size > 0;
  }

  bool is_cell_root( node const& n ) const
  {
    return _mapping_storage->mappings[this->node_to_index( n )] != 0;
  }

  void clear_mapping()
  {
    _mapping_storage->mappings.clear();
    _mapping_storage->mappings.resize( this->size(), 0 );
    _mapping_storage->mapping_size = 0;
  }

  uint32_t num_carry_lut_cells() const {
    return _mapping_storage->carry_lut_mapping_size;
  }

  uint32_t num_carry_cells() const 
  {
    return _mapping_storage->carry_mapping_size;
  }

  uint32_t num_cells() const
  {
    return _mapping_storage->mapping_size;
  }

  bool is_carry( node const& n) const {
    if (_mapping_storage->carry_mappings[this->node_to_index(n)] != 0) {
      return true;
    } else { return false; }
  }
  bool is_carry_lut( node const& n) const {
    if (n ==0) return false;
    return _mapping_storage->carry_lut_mapped[this->node_to_index(n)];
  }
  bool is_carry_lut_complemented( node const& n) const {
    return _mapping_storage->carry_lut_mapping_complemented[this->node_to_index(n)];
  }

  node carry_driver( node const& n ) const {
    if (_mapping_storage->carry_mappings[this->node_to_index(n)] != 0) {
      return _mapping_storage->carry_mappings[this->node_to_index(n)];
    } else { return 0; }
  }

  node carry_luta ( node const& n ) const {
    if (_mapping_storage->carry_luta_mappings[this->node_to_index(n)] != 0) {
      return _mapping_storage->carry_luta_mappings[this->node_to_index(n)];
    } else { return 0; }
  }

  node carry_lutb ( node const& n ) const {
    if (_mapping_storage->carry_lutb_mappings[this->node_to_index(n)] != 0) {
      return _mapping_storage->carry_lutb_mappings[this->node_to_index(n)];
    } else { return 0; }
  }

  // carry index holds the driver for this carry node
  void add_to_carry_mapping( node const& n, uint32_t carry_index )
  {
    assert (this->node_to_index( n ) < _mapping_storage->carry_mappings.size());
    assert (_mapping_storage->carry_mappings[this->node_to_index( n )] == 0);

    // carry node index
    _mapping_storage->carry_mappings[this->node_to_index( n )] = carry_index;
    _mapping_storage->carry_mapping_size++;
  }

  // carry index holds the driver for this carry node
  void add_to_carry_LUT_mapping( node const& n, uint32_t carry_LUTa_index, uint32_t carry_LUTb_index, \
      bool carry_LUTa_complemented, bool carry_LUTb_complemented )
  {
    assert (this->node_to_index( n ) < _mapping_storage->carry_luta_mappings.size());
    assert (this->node_to_index( n ) < _mapping_storage->carry_lutb_mappings.size());
    assert (_mapping_storage->carry_luta_mappings[this->node_to_index( n )] == 0);
    assert (_mapping_storage->carry_lutb_mappings[this->node_to_index( n )] == 0);

    // carry node index
    _mapping_storage->carry_luta_mappings[this->node_to_index( n )] = carry_LUTa_index;
    _mapping_storage->carry_lutb_mappings[this->node_to_index( n )] = carry_LUTb_index;
    _mapping_storage->carry_lut_mapped[this->node_to_index( carry_LUTa_index  )] = true;
    _mapping_storage->carry_lut_mapped[this->node_to_index( carry_LUTb_index  )] = true;
    _mapping_storage->carry_lut_mapping_complemented[this->node_to_index( carry_LUTa_index  )] = carry_LUTa_complemented;
    _mapping_storage->carry_lut_mapping_complemented[this->node_to_index( carry_LUTb_index  )] = carry_LUTb_complemented;
    /*std::cout << "carry node " << this->node_to_index( n ) << ": ";
    std::cout << _mapping_storage->carry_lut_mapped[this->node_to_index( carry_LUTa_index)] << " ";
    std::cout << _mapping_storage->carry_lut_mapped[this->node_to_index( carry_LUTb_index)] << " ";
    std::cout << _mapping_storage->carry_lut_mapping_complemented[this->node_to_index( carry_LUTa_index)] << " ";
    std::cout << _mapping_storage->carry_lut_mapping_complemented[this->node_to_index( carry_LUTb_index)] << " ";
    std::cout << "\n";*/
    _mapping_storage->carry_lut_mapping_size+=2;
  }

  template<typename LeavesIterator>
  void add_to_mapping( node const& n, LeavesIterator begin, LeavesIterator end )
  {

    auto& mindex = _mapping_storage->mappings[this->node_to_index( n )];

    /* increase mapping size? */
    if ( mindex == 0 )
    {
      _mapping_storage->mapping_size++;
    }

    /* set starting index of leafs */
    mindex = static_cast<uint32_t>( _mapping_storage->mappings.size() );

    /* insert number of leafs */
    _mapping_storage->mappings.push_back( static_cast<uint32_t>( std::distance( begin, end ) ) );

    /* insert leaf indexes */
    while ( begin != end )
    {
      _mapping_storage->mappings.push_back( this->node_to_index( *begin++ ) );
    }
  }

  void remove_from_mapping( node const& n )
  {
    auto& mindex = _mapping_storage->mappings[this->node_to_index( n )];

    if ( mindex != 0 )
    {
      _mapping_storage->mapping_size--;
    }

    _mapping_storage->mappings[this->node_to_index( n )] = 0;
  }

  template<bool enabled = StoreFunction, typename = std::enable_if_t<std::is_same_v<Ntk, Ntk> && enabled>>
  kitty::dynamic_truth_table cell_function( node const& n ) const
  {
    return _mapping_storage->cache[_mapping_storage->functions[this->node_to_index( n )]];
  }

  template<bool enabled = StoreFunction, typename = std::enable_if_t<std::is_same_v<Ntk, Ntk> && enabled>>
  void set_cell_function( node const& n, kitty::dynamic_truth_table const& function )
  {
    _mapping_storage->functions[this->node_to_index( n )] = _mapping_storage->cache.insert( function );
  }

  template<typename Fn>
  void foreach_cell_fanin( node const& n, Fn&& fn ) const
  {
    auto it = _mapping_storage->mappings.begin() + _mapping_storage->mappings[this->node_to_index( n )];
    const auto size = *it++;
    using IteratorType = decltype( it );
    detail::foreach_element_transform<IteratorType, typename Ntk::node>( it, it + size,
                                                                         [&]( auto i ) { return this->index_to_node( i ); }, fn );
  }

private:
  std::shared_ptr<detail::mapping_view_storage<StoreFunction>> _mapping_storage;
};

template<class T>
mapping_view(T const&) -> mapping_view<T>;

} // namespace mockturtle
