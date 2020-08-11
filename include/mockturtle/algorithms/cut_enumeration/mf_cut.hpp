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
  \file mf_cut.hpp
  \brief Cut enumeration for MF mapping (see giaMf.c)

  \author Mathias Soeken
*/

#pragma once

#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <iostream>

#include "../cut_enumeration.hpp"

namespace mockturtle
{

/*! \brief Cut implementation based on ABC's giaMf.c

  See <a href="https://github.com/berkeley-abc/abc/blob/master/src/aig/gia/giaMf.c">giaMf.c</a> in ABC's repository.
*/
struct cut_enumeration_mf_cut
{
  uint32_t delay{0};
  float flow{0};
  float cost{0};
  uint32_t carry{0};
};

template<bool ComputeTruth>
bool operator<( cut_type<ComputeTruth, cut_enumeration_mf_cut> const& c1, cut_type<ComputeTruth, cut_enumeration_mf_cut> const& c2 )
{
  constexpr auto eps{0.005f};
  if ( c1->data.flow < c2->data.flow - eps )
    return true;
  if ( c1->data.flow > c2->data.flow + eps )
    return false;
  if ( c1->data.delay < c2->data.delay )
    return true;
  if ( c1->data.delay > c2->data.delay )
    return false;
  return c1.size() < c2.size();
  //if ( c1->data.carry == 0 && c2->data.carry > 0) {
}

template<typename Cut, typename Ntk>
uint32_t count_path_to_node ( Cut& cut, Ntk const& ntk, uint32_t index, uint32_t source_index,  uint32_t dest_index) {

  auto const& n = ntk.index_to_node(source_index);

  if (source_index == dest_index) {
    return 1;
  }

  for ( auto leaf : cut )
    if (source_index == leaf) return 0;

  uint32_t total = 0;

  ntk.foreach_fanin( n, [&] ( auto const& f ){
    auto const leaf = ntk.get_node(f);
    auto const leaf_index = ntk.node_to_index(leaf);
    if (!ntk.is_constant(leaf))
      total += count_path_to_node(cut, ntk, index, leaf_index, dest_index);
  });
  return total;
}


template<>
struct cut_enumeration_update_cut<cut_enumeration_mf_cut>
{

  template<typename Cut, typename NetworkCuts, typename Ntk>
  static void apply( Cut& cut, NetworkCuts const& cuts, Ntk const& ntk, node<Ntk> const& n )
  {
    uint32_t delay{0};
    //uint32_t single_path{0};
    //uint32_t crit_path{0};
    float flow = cut->data.cost = cut.size() < 2 ? 0.0f : 1.0f;

    //auto index = ntk.node_to_index(n);

    for ( auto leaf : cut )
    {
      const auto& best_leaf_cut = cuts.cuts( leaf )[0];
      delay = std::max( delay, best_leaf_cut->data.delay );
      //if (delay == best_leaf_cut->data.delay) {
        //crit_path++;
        //if (count_path_to_node<Cut, Ntk>( cut, ntk, index, index, leaf ) == 1)
        //  single_path++;
      //}
      flow += best_leaf_cut->data.flow;
       
    }
    //cut->data.carry = crit_path - single_path;
    cut->data.delay = 1 + delay;
    cut->data.flow = flow / ntk.fanout_size( n );
  }
};

template<int MaxLeaves>
std::ostream& operator<<( std::ostream& os, cut<MaxLeaves, cut_data<false, cut_enumeration_mf_cut>> const& c )
{
  os << "{ ";
  std::copy( c.begin(), c.end(), std::ostream_iterator<uint32_t>( os, " " ) );
  os << "}, D = " << std::setw( 3 ) << c->data.delay << " A = " << c->data.flow;
  return os;
}

} // namespace mockturtle
