// Run MIG-based mapping on carry chain.

#include <string>
#include <iostream>

#include <dirent.h>

#include <mockturtle/views/mapping_view.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/carry_depth_view.hpp>
#include <mockturtle/io/verilog_reader.hpp>
#include <mockturtle/io/write_blif.hpp>
#include <mockturtle/io/write_verilog.hpp>
#include <mockturtle/algorithms/lut_mapping.hpp>
#include <mockturtle/algorithms/carry_lut_mapping.hpp>
#include <mockturtle/algorithms/collapse_mapped.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/algorithms/equivalence_checking.hpp>
#include <mockturtle/algorithms/miter.hpp>

#include <kitty/kitty.hpp>
#include <kitty/isop.hpp>

std::string getFileName(const std::string& s) {

  char dirLevel = '/';
  char period = '.';

  std::string sFile;

  size_t i = s.rfind(dirLevel, s.length());
  sFile = s.substr(i+1, s.length() - i);
  i = sFile.rfind(period);
  sFile = sFile.substr(0,i);
  i = sFile.rfind(period);
  sFile = sFile.substr(0,i);
  i = sFile.rfind(period);
  if (i != std::string::npos) {
     sFile = sFile.substr(i+1, s.length() - i);
  }
  return(sFile);
}

using namespace mockturtle;

int main (int argc, char *argv[]){

  if (argc < 2) return 0;

  //std::cout << "benchmark,nMIG,dMIG,n6LUT,d6LUT\n";

  mig_network mig;
  mig_network mig_original;
  lorina::diagnostic_engine diag;
  lorina::return_code result;

  //std::cout << argv[1] << ",";
  
  /////////////////////////////
  // Read Verilog into a MIG network.
  /////////////////////////////
  result = lorina::read_verilog(argv[1], verilog_reader(mig) ,&diag);
  if (result != lorina::return_code::success) {
    std::cout << "Parsing Error.\n"; 
  }
  result = lorina::read_verilog(argv[1], verilog_reader(mig_original) ,&diag);
  if (result != lorina::return_code::success) {
    std::cout << "Parsing Error.\n"; 
  }

  depth_view depth_mig { mig }; 
  //std::cout << mig.num_gates() << "," << depth_mig.depth() << ",";
  
  ///////////////////////////// 
  // Map MIG to LUT and carry
  /////////////////////////////
  mapping_view <mig_network, true> carry_mapped_mig { mig };
  carry_lut_mapping <mapping_view<mig_network,true>,true> (carry_mapped_mig);  
  const auto klut_carry_opt = collapse_mapped_network<klut_network>( carry_mapped_mig );
  if (klut_carry_opt == std::nullopt) {
    std::cout << "LUT and carry does not have mapping\n";
    return 0;
  }
  auto const& klut_carry = *klut_carry_opt;
  carry_depth_view depth_klut_carry { klut_carry }; 
  std::cout << "lut mapping finished\n";

  /////////////////////////////
  // Map MIG to LUT only
  /////////////////////////////
  mapping_view <mig_network, true> mapped_mig { mig };
  lut_mapping <mapping_view<mig_network,true>,true> (mapped_mig); 
  const auto klut_opt = collapse_mapped_network<klut_network>( mapped_mig );
  if (klut_opt == std::nullopt) {
    std::cout << "Does not have mapping\n";
    return 0;
  }
  auto const& klut = *klut_opt;
  depth_view depth_klut { klut }; 
  
  /////////////////////////////
  // Check circuit equivalence.
  /////////////////////////////
  //const auto miter_circuit = miter<mig_network, mig_network, mig_network> (mig, mig_original);
  //const auto miter_circuit = miter<klut_network, klut_network, mig_network> (klut_carry, mig_original);
  /*const auto miter_circuit = miter<klut_network, klut_network, klut_network> (klut_carry, klut);
  if ( miter_circuit == std::nullopt) std::cout << "Input/output numbers do not match.\n";
  const auto result_of_equivalence = equivalence_checking ( *miter_circuit );
  if ( !result_of_equivalence || ! *result_of_equivalence ) {
    std::cout << "Networks are different.\n";
    return 0;
  }*/

  //std::cout << klut.num_gates() << "," << depth_klut.depth() << ",";
  //std::cout << klut_carry.num_gates() << "," << float(depth_klut_carry.depth()/7.0) << "\n";
  write_blif(klut_carry, "blif/" + getFileName(argv[1]) + "_carry.blif");  
  write_blif(klut, "blif/" + getFileName(argv[1]) + ".blif");  

  std::ofstream os( ("blif/" + getFileName(argv[1]) + "_carry.v").c_str(), std::ofstream::out );
  write_verilog(mig, os);
  std::ofstream os1( ("blif/" + getFileName(argv[1]) + ".v").c_str(), std::ofstream::out );
  write_verilog(mig_original, os1);

  mig.clear_values();

  return 0;

}

