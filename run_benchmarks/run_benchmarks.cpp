// Run MIG-based mapping on carry chain.

#include <string>
#include <iostream>

#include <dirent.h>


#include <mockturtle/io/verilog_reader.hpp>
#include <mockturtle/io/write_blif.hpp>
#include <mockturtle/algorithms/carry_lut_mapping.hpp>
#include <mockturtle/algorithms/collapse_mapped.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/carry_depth_view.hpp>
#include <mockturtle/views/mapping_view.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/io/write_dot.hpp>

using namespace mockturtle;

std::string getFileName(const std::string& s) {
  std::string s_rtn;
  s_rtn = s.substr(s.find_last_of('/')+1);
  //s_rtn.erase(0,s_rtn.find('.')+1); 
  s_rtn.erase(s_rtn.find('.')); 
  return s_rtn;
}

int main (int argc, char *argv[]){
  mig_network mig;
  lorina::diagnostic_engine diag;
  lorina::return_code result;
  
  /////////////////////////////
  // Read Verilog into a MIG network.
  /////////////////////////////
  result = lorina::read_verilog(argv[1], verilog_reader(mig) ,&diag);
  if (result != lorina::return_code::success) {
    std::cout << "Parsing Error.\n"; 
  }
  depth_view depth_mig { mig }; 

  //////////////////////////////////////////////////////
  // Map MIG to LUT and carry
  //////////////////////////////////////////////////////
  // argv[1]        benchmark name
  // argv[2](true)  do carry mapping?
  // argv[3](true)  use xilinx arch?
  // argv[4](300)   number of max carry mapping
  // argv[5](4)     cost function used 

  std::string outputName; 

  carry_lut_mapping_params mapping_params;

  // Set carry LUT mapping parameters
  if (argc > 3 && argv[3] == std::string("delay"))
    mapping_params.delay = true;
  if (argc > 4 && argv[4] == std::string("xilinx"))
    mapping_params.xilinx_arch = true;
  if (argc > 2 && argv[2] == std::string("baseline")) {
    mapping_params.carry_mapping = false;
    outputName = getFileName(argv[1]);
  } else {
    mapping_params.carry_mapping = true;
    outputName = getFileName(argv[1]) + "_carry";
    if (argc > 5)
      mapping_params.max_rounds_carry = std::stoi(argv[5]); 
    if (argc > 6)
      mapping_params.cost = std::stoi(argv[6]);
  }

  mapping_view <mig_network, true> carry_mapped_mig { mig };
  carry_lut_mapping <mapping_view<mig_network,true>,true> (carry_mapped_mig, mapping_params);  
  auto klut_carry_opt = collapse_mapped_network<klut_network>( carry_mapped_mig );
  if (klut_carry_opt == std::nullopt) {
    std::cout << "LUT and carry does not have mapping\n";
    return 0;
  }
  auto& klut_carry = *klut_carry_opt;
  carry_depth_view depth_klut_carry { klut_carry }; 

  klut_carry.foreach_node( [&]( auto const& n ) {
    if ( klut_carry.is_constant( n ) || klut_carry.is_pi( n ) || klut_carry.is_carry(n)  )
      return; /* continue */
    if (klut_carry.is_carry_LUT(n)) {
      klut_carry.increase_carry_LUT_savings();
      return;
    }
  });
 
  // Write blif output. Last line will be a comment on what run it was 
  write_blif(klut_carry, "blif/" + outputName+".blif", true/*carry mapping*/, mapping_params.xilinx_arch );
  std::ofstream blifOut ("blif/" + outputName+".blif", std::ofstream::out | std::ofstream::app );
  blifOut << "# " << outputName << " baseline=" << !mapping_params.carry_mapping
          << " xilinx=" << mapping_params.xilinx_arch << " maxCarryRounds="
          << mapping_params.max_rounds_carry  << " cost=" << mapping_params.cost << "\n\n"; 
  blifOut.close();

  std::cout << "Results for " << outputName  << ",";
  std::cout << mig.num_gates() << "," << depth_mig.depth();
  //std::cout << ",=" << klut_carry.num_gates() << "+" << klut_carry.num_carry() << "/2," << float(depth_klut_carry.depth()/LUT_DELAY);
  std::cout << "," << klut_carry.num_gates() + klut_carry.num_carry() - klut_carry.num_carry_LUT(); 
  std::cout << "," << klut_carry.num_carry();
  std::cout << "," << klut_carry.num_carry_LUT() << "," << float(depth_klut_carry.depth()/LUT_DELAY);
  std::cout << "\n";

  //depth_mig.print_levels();
  //depth_klut_carry.print_levels();
  //depth_klut_carry.print_num_paths_at_max_level();

  // Write blif for CEC. Doesn't have to happen every run.
  write_blif(klut_carry, "blifcec/" + outputName +".cec.blif", false/*carry mapping*/, mapping_params.xilinx_arch );

  mig.clear_values();

  return 0;

}

