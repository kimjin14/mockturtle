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
//#include <mockturtle/algorithms/carry_chain_mapping.hpp>
#include <mockturtle/algorithms/collapse_mapped.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/networks/xmg.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/algorithms/equivalence_checking.hpp>
#include <mockturtle/algorithms/miter.hpp>
#include <mockturtle/io/write_dot.hpp>

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

  bool runbaseline = false;
  if (argc > 2 && argv[2] == std::string("baseline")) {
    std::cout << "Baseline.\n";
    runbaseline = true;
  }

  mig_network mig;
  mig_network mig_original;
  lorina::diagnostic_engine diag;
  lorina::return_code result;
  
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
  //xmg_network xmg;
  //result = lorina::read_verilog(argv[1], verilog_reader(xmg) ,&diag);
  //if (result != lorina::return_code::success) {
  //  std::cout << "Parsing Error.\n"; 
  //}
  //depth_view depth_xmg { xmg };

  //std::cout << "depth of xmg is " << depth_xmg.depth() << "\n";

  
  /////////////////////////////
  // Map MIG to LUT only
  /////////////////////////////
  if (runbaseline) {
    mapping_view <mig_network, true> mapped_mig { mig };
    lut_mapping <mapping_view<mig_network,true>,true>(mapped_mig); 
    const auto klut_opt = collapse_mapped_network<klut_network>( mapped_mig );
    if (klut_opt == std::nullopt) {
      std::cout << "Does not have mapping\n";
      return 0;
    }
    auto const& klut = *klut_opt;
    depth_view depth_klut ( klut ); 

    //write_dot(klut, "output_lut.dot");
    write_blif(klut, "blif/" + getFileName(argv[1]) + ".blif");  

    std::cout << "Results for " << argv[1] << ",";
    std::cout << mig.num_gates() << "," << depth_mig.depth();
    std::cout << "," << klut.num_gates() << "," << depth_klut.depth();
    std::cout << "\n";
  }
  ///////////////////////////// 
  // Map MIG to LUT and carry
  /////////////////////////////
  else {
    std::cout << "Creating mapping view\n";
    mapping_view <mig_network, true> carry_mapped_mig { mig };
    carry_lut_mapping <mapping_view<mig_network,true>,true> (carry_mapped_mig);  
    const auto klut_carry_opt = collapse_mapped_network<klut_network>( carry_mapped_mig );
    if (klut_carry_opt == std::nullopt) {
      std::cout << "LUT and carry does not have mapping\n";
      return 0;
    }
    auto const& klut_carry = *klut_carry_opt;

    carry_depth_view depth_klut_carry { klut_carry }; 

    //write_dot(klut_carry, "output_carry_lut.dot");
    write_blif(klut_carry, "blif/" + getFileName(argv[1]) + "_carry.blif");  

    std::cout << "Results for " << argv[1] << ",";
    std::cout << mig.num_gates() << "," << depth_mig.depth();
    std::cout << ",=" << klut_carry.num_gates() << "+" << klut_carry.num_carry() << "/2," << float(depth_klut_carry.depth()/LUT_DELAY);
    std::cout << "\n";
  }

  /////////////////////////////
  // Print mapping results.
  /////////////////////////////
  //std::cout << "Results for " << argv[1] << ",";
  //std::cout << mig.num_gates() << "," << depth_mig.depth();
  //std::cout << "," << klut_carry.num_gates() << "(" << klut_carry.num_carry() << ")," << float(depth_klut_carry.depth()/LUT_DELAY);
  //std::cout << "\n";

  /////////////////////////////
  // Write to Verilog for backup
  /////////////////////////////
  //std::ofstream os( ("blif/" + getFileName(argv[1]) + "_carry.v").c_str(), std::ofstream::out );
  //write_verilog(mig, os);
  //std::ofstream os1( ("blif/" + getFileName(argv[1]) + ".v").c_str(), std::ofstream::out );
  //write_verilog(mig_original, os1);

  mig.clear_values();

  return 0;

}

