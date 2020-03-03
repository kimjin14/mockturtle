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
  \file write_blif.hpp
  \brief Write networks to BLIF format

  \author Heinz Riener
*/

#pragma once

#include <fstream>
#include <iostream>
#include <string>

#include <fmt/format.h>
#include <kitty/operations.hpp>
#include <kitty/isop.hpp>

#include "../traits.hpp"
#include "../views/topo_view.hpp"

#define CARRY_MAPPING false 

namespace mockturtle
{

/*! \brief Writes network in BLIF format into output stream
 *
 * An overloaded variant exists that writes the network into a file.
 *
 * **Required network functions:**
 * - `fanin_size`
 * - `foreach_fanin`
 * - `foreach_pi`
 * - `foreach_po`
 * - `get_node`
 * - `is_constant`
 * - `is_pi`
 * - `node_function`
 * - `node_to_index`
 * - `num_pis`
 * - `num_pos`
 *
 * \param ntk Network
 * \param os Output stream
 */
template<class Ntk>
void write_blif( Ntk const& ntk, std::ostream& os  )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_fanin_size_v<Ntk>, "Ntk does not implement the fanin_size method" );
  static_assert( has_foreach_fanin_v<Ntk>, "Ntk does not implement the foreach_fanin method" );
  static_assert( has_foreach_pi_v<Ntk>, "Ntk does not implement the foreach_pi method" );
  static_assert( has_foreach_po_v<Ntk>, "Ntk does not implement the foreach_po method" );
  static_assert( has_is_constant_v<Ntk>, "Ntk does not implement the is_constant method" );
  static_assert( has_is_pi_v<Ntk>, "Ntk does not implement the is_pi method" );
  static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
  static_assert( has_num_pis_v<Ntk>, "Ntk does not implement the num_pis method" );
  static_assert( has_num_pos_v<Ntk>, "Ntk does not implement the num_pos method" );
  static_assert( has_node_to_index_v<Ntk>, "Ntk does not implement the node_to_index method" );
  static_assert( has_node_function_v<Ntk>, "Ntk does not implement the node_function method" );

  topo_view topo_ntk{ntk};

  // This file contains informations on which nodes are carry chains
  // and which blif outputs are related to those carry chains

  os << ".model netlist\n";

  if ( topo_ntk.num_pis() > 0u )
  {
    os << ".inputs ";
    uint32_t count = 0;
    topo_ntk.foreach_pi( [&]( auto const& n ) {
      os << fmt::format( "n{} ", topo_ntk.node_to_index( n ));
      if (count > 100) {
        os << fmt::format( "\\\n" );
        count = 0;
      } else
        count++;
    } );
    os << "\n";
  }

  if ( topo_ntk.num_pos() > 0u )
  {
    os << ".outputs ";
    uint32_t count = 0;
    topo_ntk.foreach_po( [&]( auto const& n, auto index ) {
        os << fmt::format( "po{} ", index );
        if (count > 100) {
          os << fmt::format( "\\\n" );
          count = 0;
        } else
          count++;
      } );
    os << "\n";
  }

  os << ".names n0\n";
  os << "0\n";

  os << ".names n1\n";
  os << "1\n";

  // Place all nodes into a vector first
  // Number of list can change depending on how many carry chains
  std::vector< std::vector<uint32_t> >carry_chains;

  std::vector<uint32_t> adder_cout;
  std::vector<uint32_t> adder_a;
  std::vector<uint32_t> adder_b;
  std::vector<uint32_t> adder_cin;

  uint32_t separate_chain = 0;

  uint32_t next_node = topo_ntk.size();
  topo_ntk.foreach_node( [&]( auto const& n ) {

    if ( topo_ntk.is_constant( n ) || topo_ntk.is_pi( n ) )
      return; /* continue */

    // In case of carry mapped LUT
    if (CARRY_MAPPING == true && topo_ntk.is_carry(n) ) {
  
      auto const func = topo_ntk.node_function( n );
      auto list_of_cubes = isop(func);     
      uint32_t count = 0;

      uint32_t should_be_constant_0 = 0;
      uint32_t should_be_constant_1 = 0;

      auto const child = topo_ntk.node_to_index(topo_ntk.get_children(n, topo_ntk.fanin_size(n)-1) );
      bool inserted = false; 
        for (auto& carry_chain: carry_chains) {
          if (!carry_chain.empty() && carry_chain[carry_chain.size()-1] == child) {
            carry_chain.push_back(n);
            inserted = true;
          }
        }
      if ( !inserted ) {
          std::vector<uint32_t> new_chain;
          new_chain.push_back(n);
          carry_chains.push_back(new_chain);
      } 
    } else {

      os << fmt::format( ".names " );
      topo_ntk.foreach_fanin( n, [&]( auto const& c ) {
        os << fmt::format( "n{} ", topo_ntk.node_to_index( c ) );
      });
      
      os << fmt::format( " n{}\n", topo_ntk.node_to_index( n ) );
    
      auto const func = topo_ntk.node_function( n );
      for ( const auto& cube : isop( func ) )
      {
        cube.print( topo_ntk.fanin_size( n ), os );
        os << " 1\n";
      }
    }
  });

  for (auto carry_chain: carry_chains) {
    separate_chain = 0;
    next_node = topo_ntk.size();
    for (auto n: carry_chain) {

      if (separate_chain == 10) {
        os << fmt::format(".names n{} n{}\n", topo_ntk.node_to_index(topo_ntk.get_children(n, topo_ntk.fanin_size(n)-1) ), next_node);
        os << fmt::format("1 1\n\n", topo_ntk.node_to_index(topo_ntk.get_children(n, topo_ntk.fanin_size(n)-1) ), next_node);
        separate_chain = 0;
      }

      os << (".subckt adder_lut ");

      uint32_t input_count = 0;
      topo_ntk.foreach_fanin( n, [&]( auto const& c ) {
        if (input_count != topo_ntk.fanin_size(n)-1) {
          os << fmt::format( "in{}=n{} ", input_count, topo_ntk.node_to_index( c ) );
          input_count++;
        }
      });
      for (; input_count < 5; input_count++)
        os << fmt::format( "in{}=unconn ", input_count );
      
      assert(topo_ntk.fanin_size(n) <= 6); 
      if (topo_ntk.is_pi(topo_ntk.get_children(n, topo_ntk.fanin_size(n)-1))) {
        os << fmt::format("cin=n{} ", topo_ntk.node_to_index(topo_ntk.get_children(n, topo_ntk.fanin_size(n)-1) ) );
      } else if (separate_chain == 0) { 
        os << fmt::format("cin=n{} ", next_node++);
      } else {
        os << fmt::format("cin=c{} ", topo_ntk.node_to_index(topo_ntk.get_children(n, topo_ntk.fanin_size(n)-1) ) );
      }

      os << fmt::format("cout=c{} ", topo_ntk.node_to_index( n ) );
      os << fmt::format("sumout=n{}\n", topo_ntk.node_to_index( n ) );

      separate_chain++;
    }
    os << fmt::format("\n");
  }
  os << fmt::format("\n");

  if ( topo_ntk.num_pos() > 0u )
  {
    topo_ntk.foreach_po( [&]( auto const& n, auto index ) {
          os << fmt::format( ".names n{} po{}\n1 1\n", topo_ntk.node_to_index( n ), index );
      } );
  }
  os << ".end\n";

  if (CARRY_MAPPING) {
    os << "\n.model adder_lut\n";
    os << ".inputs in0 in1 in2 in3 in4 cin\n";
    os << ".outputs sumout cout\n";
    os << ".blackbox\n";
    os << ".end\n";
  }
  os << std::flush;
}

/*! \brief Writes network in BENCH format into a file
 *
 * **Required network functions:**
 * - `fanin_size`
 * - `foreach_fanin`
 * - `foreach_pi`
 * - `foreach_po`
 * - `get_node`
 * - `is_constant`
 * - `is_pi`
 * - `node_function`
 * - `node_to_index`
 * - `num_pis`
 * - `num_pos`
 *
 * \param ntk Network
 * \param filename Filename
 */
template<class Ntk>
void write_blif( Ntk const& ntk, std::string const& filename )
{
  std::ofstream os( filename.c_str(), std::ofstream::out );
  write_blif( ntk, os );
  os.close();
}

} /* namespace mockturtle */