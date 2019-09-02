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
void write_blif( Ntk const& ntk, std::ostream& os, std::ostream& os_log  )
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
  os_log << "Carry chain mapping info.\n";

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

  uint32_t next_node = topo_ntk.size();

  std::vector<uint32_t> adder_cout;
  std::vector<uint32_t> adder_a;
  std::vector<uint32_t> adder_b;
  std::vector<uint32_t> adder_cin;

  topo_ntk.foreach_node( [&]( auto const& n ) {
    os_log << n << ":\n";    
    if ( topo_ntk.is_constant( n ) || topo_ntk.is_pi( n ) )
      return; /* continue */

    // In case of carry mapped LUT
    if (topo_ntk.is_carry(n) ) {
  
      auto const func = topo_ntk.node_function( n );
      auto list_of_cubes = isop(func);     
      uint32_t count = 0;

      uint32_t should_be_constant_0 = 0;
      uint32_t should_be_constant_1 = 0;

      // The inputs to the carry node should be separated into
      // 2 separate LUT mapping, which doesn't include the carry-in
      // from the previous carry node
      for ( const auto& cube : list_of_cubes ) {
        if(!cube.get_mask(topo_ntk.fanin_size(n) - 1)) {
          if (cube.num_literals() == 0) should_be_constant_1 ++;
          should_be_constant_0++;
        } else if (cube.get_bit(topo_ntk.fanin_size(n) - 1) == 0) {
          if (cube.num_literals() == 1) should_be_constant_1 ++;
          should_be_constant_0++;
        }
      }

      os << fmt::format( ".names " );
      if (should_be_constant_0 > 0 && should_be_constant_1 == 0) {
        count = 0;
        topo_ntk.foreach_fanin( n, [&]( auto const& c ) {
          if (count != topo_ntk.fanin_size(n)-1) // don't include carry-in
            os << fmt::format( "n{} ", topo_ntk.node_to_index( c ) );
          count++;
        });
        os << fmt::format( " n{}\n", next_node++ );
        for ( const auto& cube : list_of_cubes ) {
          if (!cube.get_mask(topo_ntk.fanin_size(n) - 1) || cube.get_bit(topo_ntk.fanin_size(n) - 1) == 0) {
            cube.print( topo_ntk.fanin_size( n )-1, os );
            os << " 1\n";
          }
        }
      } else {
        if (should_be_constant_1 != 0) {
          os << fmt::format( "n1 n{}\n", next_node++ );
          os << "1 1\n";
        } else {
          os << fmt::format( "n0 n{}\n", next_node++ );
          os << "1 1\n";
        }
      }
      should_be_constant_0 = 0;
      should_be_constant_1 = 0;
      for ( const auto& cube : list_of_cubes ) {
        if(!cube.get_mask(topo_ntk.fanin_size(n) - 1)) {
          if (cube.num_literals() == 0) should_be_constant_1 ++;
          should_be_constant_0++;
        } else if (cube.get_bit(topo_ntk.fanin_size(n) - 1) == 1) {
          if (cube.num_literals() == 1) should_be_constant_1 ++;
          should_be_constant_0++;
        }
      }
      os << fmt::format( ".names " );
      if (should_be_constant_0 > 0 && should_be_constant_1 == 0) {
        count = 0;
        topo_ntk.foreach_fanin( n, [&]( auto const& c ) {
          if (count != topo_ntk.fanin_size(n)-1) // don't include carry-in
            os << fmt::format( "n{} ", topo_ntk.node_to_index( c ) );
          count++;
        });
        os << fmt::format( " n{}\n", next_node++ );
        for ( const auto& cube : list_of_cubes ) {
          if (!cube.get_mask(topo_ntk.fanin_size(n) - 1) || cube.get_bit(topo_ntk.fanin_size(n) - 1) == 1) {
            cube.print( topo_ntk.fanin_size( n )-1, os );
            os << " 1\n";
          }
        }
      } else {
        if (should_be_constant_1 != 0) {
          os << fmt::format( "n1 n{}\n", next_node++ );
          os << "1 1\n";
        } else {
          os << fmt::format( "n0 n{}\n", next_node++ );
          os << "1 1\n";
        }
      }

      auto carry_child = topo_ntk.get_children(n, topo_ntk.fanin_size(n)-1);
      adder_a.push_back(next_node-2);
      adder_b.push_back(next_node-1);
      adder_cin.push_back(topo_ntk.node_to_index(carry_child));
      adder_cout.push_back(topo_ntk.node_to_index(n));

    } else {

      os << fmt::format( ".names " );
      topo_ntk.foreach_fanin( n, [&]( auto const& c ) {
        os << fmt::format( "n{} ", topo_ntk.node_to_index( c ) );
      });
      if (topo_ntk.is_carry(n)) 
        os << fmt::format( " n{}\n", topo_ntk.node_to_index( n ) );
      else
        os << fmt::format( " n{}\n", topo_ntk.node_to_index( n ) );
    
      auto const func = topo_ntk.node_function( n );
      for ( const auto& cube : isop( func ) )
      {
        cube.print( topo_ntk.fanin_size( n ), os );
        os << " 1\n";
      }
    }
  });

  for (uint32_t i = 0; i < adder_cout.size(); i++) {
    // Print mig node
    if (true) {
      os << fmt::format( ".subckt adder " );
      os << fmt::format( "a=n{} ", adder_a[i]);
      os << fmt::format( "b=n{} ", adder_b[i]);
      os << fmt::format( "cin=c{} ", adder_cin[i]);
      os << fmt::format( "sumout=n{} ",adder_cout[i]);
      os << fmt::format( "cout=c{}\n", adder_cout[i]);
    } else { 
      os << fmt::format( ".names " );
      os << fmt::format( "n{} ", adder_a[i]);
      os << fmt::format( "n{} ", adder_b[i]);
      os << fmt::format( "n{} ", adder_cin[i]);
      os << fmt::format( " n{}\n", adder_cout[i]);
      os << "-11 1\n";
      os << "11- 1\n";
      os << "1-1 1\n";
    }
    os_log << "\t" << "carry LUT nodes ";    
    os_log << "n" << adder_a[i] << " ";
    os_log << "n" << adder_b[i] << "\n";

  }

  if ( topo_ntk.num_pos() > 0u )
  {
    topo_ntk.foreach_po( [&]( auto const& n, auto index ) {
        //if (topo_ntk.is_carry(n) && topo_ntk.is_complemented(n))
        //  os << fmt::format( ".names n{} po{}\n0 1\n", topo_ntk.node_to_index( n ), index );
        //else
          os << fmt::format( ".names n{} po{}\n1 1\n", topo_ntk.node_to_index( n ), index );
      } );
  }
  os << ".end\n";


  os << "\n.model adder\n";
  os << ".inputs a b cin\n";
  os << ".outputs sumout cout\n";
  os << ".blackbox\n";
  os << ".end\n\n";

  os << std::flush;
  os_log << std::flush;
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
void write_blif( Ntk const& ntk, std::string const& filename, std::string const& filename_log )
{
  std::ofstream os_log( filename_log.c_str(), std::ofstream::out );
  std::ofstream os( filename.c_str(), std::ofstream::out );
  write_blif( ntk, os, os_log );
  os.close();
  os_log.close();
}

} /* namespace mockturtle */
