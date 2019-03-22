// Run MIG-based mapping on carry chain.

#include <string>
#include <iostream>

#include <dirent.h>

#include <mockturtle/views/mapping_view.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/io/verilog_reader.hpp>
#include <mockturtle/algorithms/carry_lut_mapping.hpp>
#include <mockturtle/algorithms/collapse_mapped.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/networks/klut.hpp>

#include <kitty/kitty.hpp>

using namespace mockturtle;

int main (int argc, char *argv[]){

  if (argc < 2) return 0;

  std::cout << "benchmark,nMIG,dMIG,n6LUT,d6LUT\n";

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
  carry_lut_mapping (mapped_mig);  
  //lut_mapping <mapping_view<mig_network,true>,true> (mapped_mig); // for storing fcn
  //lut_mapping (mapped_mig);

  // Collapse mapped MIG to LUT network
  const auto klut_opt = collapse_mapped_network<klut_network>( mapped_mig );
  auto const& klut = *klut_opt;
  depth_view depth_klut { klut }; 

  std::cout << klut.num_gates() << "," << depth_klut.depth();
  //mig.clear_values();

  std::cout << "\n";

  return 0;

}
