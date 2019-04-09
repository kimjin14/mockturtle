// Run MIG-based mapping on carry chain.

#include <string>
#include <iostream>

#include <dirent.h>

#include <mockturtle/views/mapping_view.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/carry_depth_view.hpp>
#include <mockturtle/io/verilog_reader.hpp>
#include <mockturtle/io/write_verilog.hpp>
//#include <mockturtle/algorithms/lut_mapping.hpp>
#include <mockturtle/algorithms/carry_lut_mapping.hpp>
#include <mockturtle/algorithms/collapse_mapped.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/networks/klut.hpp>

#include <kitty/kitty.hpp>

using namespace mockturtle;

int main (int argc, char *argv[]){

  if (argc < 2) return 0;

  //std::cout << "benchmark,nMIG,dMIG,n6LUT,d6LUT\n";

  mig_network mig;
  lorina::diagnostic_engine diag;
  lorina::return_code result;

  //std::cout << argv[1] << ",";
  
  // Read Verilog into a MIG network.
  result = lorina::read_verilog(argv[1], verilog_reader(mig) ,&diag);
  if (result != lorina::return_code::success) {
    std::cout << "Parsing Error.\n"; 
  }

  // Use depth_view to report the depth of the MIG network.
  depth_view depth_mig { mig }; 
  //std::cout << mig.num_gates() << "," << depth_mig.depth() << ",";
 
  // Map MIG to 6LUT
  mapping_view <mig_network, true> mapped_mig { mig };
  //carry_lut_mapping (mapped_mig);  
  carry_lut_mapping <mapping_view<mig_network,true>,false> (mapped_mig);  
  //lut_mapping <mapping_view<mig_network,true>,true> (mapped_mig); // for storing fcn
  //lut_mapping (mapped_mig);

  //std::cout << mapped_mig.num_carry_cells() << "," << mapped_mig.num_carry_lut_cells() << ",";

  // Collapse mapped MIG to LUT network
  const auto klut_opt = collapse_mapped_network<klut_network>( mapped_mig );
  if (klut_opt == std::nullopt) {
    std::cout << "Does not have mapping\n";
    return 0;
  }

  auto const& klut = *klut_opt;

 /* klut.foreach_gate([&](auto n) {
    if (klut.is_carry(n)) {
      std::cout << "carry ";
    } else {
      std::cout << "regular ";
    }
    std::cout << n << ":";
    kitty::print_hex(klut.node_function(n));
    std::cout << "\n";
  });*/

  carry_depth_view depth_klut { klut }; 

  std::cout << klut.num_gates() << "," << float(depth_klut.depth()/7.0);
  mig.clear_values();

  std::cout << "\n";

  //write_verilog(klut, "test.v");

  return 0;

}

