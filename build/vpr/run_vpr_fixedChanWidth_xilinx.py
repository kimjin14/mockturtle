#!/usr/bin/env python3

from shutil import copyfile
import os
import sys

channelWidthDict = {
  "mem_ctrl":89,"usb_phy":38,"diffeq1":76,"systemcaes":81,"MUL32":90,"MAC32":94,"sasc":39,"i2c":42,"hamming":74,"div16":86,"max":66,"tv80":98,"rca_128":146,"rca_32":60,"spi":70,"revx":71,"CSA464":56,"log2":168,"systemcdes":47,"square":109,"simple_spi":38,"tern_add":66,"mult64":124,"rca_256":68,"aes_core":64,"comp":178,"pci_spoci_ctrl":65,"sqrt32":75,"rca_64":55,"des_area":79,"ss_pcm":37,
  "mem_ctrl_carry":89,"usb_phy_carry":38,"diffeq1_carry":76,"systemcaes_carry":81,"MUL32_carry":90,"MAC32_carry":94,"sasc_carry":39,"i2c_carry":42,"hamming_carry":74,"div16_carry":86,"max_carry":66,"tv80_carry":98,"rca_128_carry":146,"rca_32_carry":60,"spi_carry":70,"revx_carry":71,"CSA464_carry":56,"log2_carry":168,"systemcdes_carry":47,"square_carry":109,"simple_spi_carry":38,"tern_add_carry":66,"mult64_carry":124,"rca_256_carry":68,"aes_core_carry":64,"comp_carry":178,"pci_spoci_ctrl_carry":65,"sqrt32_carry":75,"rca_64_carry":55,"des_area_carry":79,"ss_pcm_carry":37,
  "depth_i2c":49,"area_sin":115,"area_int2float":24,"depth_bar":70,"depth_sin":118,"area_cavlc":28,"area_bar":64,"area_i2c":45,"area_dec":40,"area_voter":49,"area_priority":60,"depth_int2float":24,"area_max":75,"depth_cavlc":32,"area_ctrl":42,"area_arbiter":148,"area_log2":166,"area_square":111,"depth_ctrl":43,"depth_adder":55,"depth_dec":40,"depth_arbiter":150,"depth_max":78,"depth_log2":176,"depth_voter":56,"area_router":49,"area_mem_ctrl":123,"area_multiplier":126,"depth_multiplier":138,"depth_priority":62,"area_adder":51,"depth_square":118,"depth_mem_ctrl":145,"depth_router":41,
  "depth_i2c_carry":49,"area_sin_carry":115,"area_int2float_carry":24,"depth_bar_carry":70,"depth_sin_carry":118,"area_cavlc_carry":28,"area_bar_carry":64,"area_i2c_carry":45,"area_dec_carry":40,"area_voter_carry":49,"area_priority_carry":60,"depth_int2float_carry":24,"area_max_carry":75,"depth_cavlc_carry":32,"area_ctrl_carry":42,"area_arbiter_carry":148,"area_log2_carry":166,"area_square_carry":111,"depth_ctrl_carry":43,"depth_adder_carry":55,"depth_dec_carry":40,"depth_arbiter_carry":150,"depth_max_carry":78,"depth_log2_carry":176,"depth_voter_carry":56,"area_router_carry":49,"area_mem_ctrl_carry":123,"area_multiplier_carry":126,"depth_multiplier_carry":138,"depth_priority_carry":62,"area_adder_carry":51,"depth_square_carry":118,"depth_mem_ctrl_carry":145,"depth_router_carry":41,
  # same as area
  "preopt_adder":51,"preopt_arbiter":148,"preopt_bar":64,"preopt_cavlc":28,"preopt_ctrl":42,"preopt_dec":40,"preopt_i2c":45,"preopt_int2float":24,"preopt_log2":166,"preopt_max":75,"preopt_mem_ctrl":123,"preopt_multiplier":126,"preopt_priority":60,"preopt_router":49,"preopt_sin":115,"preopt_square":111,"preopt_voter":49,
  "preopt_adder_carry":51,"preopt_arbiter_carry":148,"preopt_bar_carry":64,"preopt_cavlc_carry":28,"preopt_ctrl_carry":42,"preopt_dec_carry":40,"preopt_i2c_carry":45,"preopt_int2float_carry":24,"preopt_log2_carry":166,"preopt_max_carry":75,"preopt_mem_ctrl_carry":123,"preopt_multiplier_carry":126,"preopt_priority_carry":60,"preopt_router_carry":49,"preopt_sin_carry":115,"preopt_square_carry":111,"preopt_voter_carry":49

#  "area_sin":117,"area_int2float":25,"area_cavlc":28,"area_bar":64,"area_i2c":45,"area_dec":40,"area_voter":52,"area_priority":58,"area_max":75,"area_ctrl":41,"area_arbiter":150,"area_log2":169,"area_square":110,"area_router":51,"area_mem_ctrl":127,"area_multiplier":126,"area_adder":49,
#  "area_sin_carry":117,"area_int2float_carry":25,"area_cavlc_carry":28,"area_bar_carry":64,"area_i2c_carry":45,"area_dec_carry":40,"area_voter_carry":52,"area_priority_carry":58,"area_max_carry":75,"area_ctrl_carry":41,"area_arbiter_carry":150,"area_log2_carry":169,"area_square_carry":110,"area_router_carry":51,"area_mem_ctrl_carry":127,"area_multiplier_carry":126,"area_adder_carry":49,
#  "depth_i2c":48,"depth_bar":69,"depth_sin":119,"depth_int2float":24,"depth_cavlc":31,"depth_ctrl":43,"depth_adder":53,"depth_dec":40,"depth_arbiter":147,"depth_max":76,"depth_log2":174,"depth_voter":59,"depth_multiplier":138,"depth_priority":62,"depth_square":117,"depth_mem_ctrl":152,"depth_router":42,
#  "depth_i2c_carry":48,"depth_bar_carry":69,"depth_sin_carry":119,"depth_int2float_carry":24,"depth_cavlc_carry":31,"depth_ctrl_carry":43,"depth_adder_carry":53,"depth_dec_carry":40,"depth_arbiter_carry":147,"depth_max_carry":76,"depth_log2_carry":174,"depth_voter_carry":59,"depth_multiplier_carry":138,"depth_priority_carry":62,"depth_square_carry":117,"depth_mem_ctrl_carry":152,"depth_router_carry":42,
#  "aes_core":62, "comp":178, "CSA464":62, "des_area":80, "diffeq1":76, "div16":86, "hamming":76, "i2c":42, "log2":162, "MAC32":90, "max":62, "mem_ctrl":88, "MUL32":94, "mult64":124, "pci_spoci_ctrl":70, "rca_128":148, "rca_256":70, "rca_32":62, "rca_64":56, "revx":70, "sasc":38, "simple_spi":40, "spi":66, "sqrt32":76, "square":110, "ss_pcm":36, "systemcaes":86, "systemcdes":50, "tern_add":70, "tv80":98, "usb_phy":38,
#  "aes_core_carry":62, "comp_carry":178, "CSA464_carry":62, "des_area_carry":80, "diffeq1_carry":76, "div16_carry":86, "hamming_carry":76, "i2c_carry":42, "log2_carry":162, "MAC32_carry":90, "max_carry":62, "mem_ctrl_carry":88, "MUL32_carry":94, "mult64_carry":124, "pci_spoci_ctrl_carry":70, "rca_128_carry":148, "rca_256_carry":70, "rca_32_carry":62, "rca_64_carry":56, "revx_carry":70, "sasc_carry":38, "simple_spi_carry":40, "spi_carry":66, "sqrt32_carry":76, "square_carry":110, "ss_pcm_carry":36, "systemcaes_carry":86, "systemcdes_carry":50, "tern_add_carry":70, "tv80_carry":98, "usb_phy_carry":38
}

blifDir= "/home/kimjin14/work/mapping/mig_mapping/mockturtle/build/"
vprDir="/home/kimjin14/work/vtr/vtr-verilog-to-routing/vpr/"
runDir="/home/kimjin14/work/mapping/mig_mapping/mockturtle/build/vpr/"
outputDir="/home/kimjin14/work/mapping/mig_mapping/mockturtle/build/vpr/output"
saveDir="/home/kimjin14/work/mapping/mig_mapping/mockturtle/build/vpr/output"

xmlFile = "arch/Xilinx_style_arch_ecin.xml"
vprOptions = "--timing_report_detail aggregated --route_chan_width"

# Get directory of what you want to run
blifDir=blifDir+sys.argv[1]+"/"

#seed_list={1, 2, 5}
seed_list={1, 5}

for seed in seed_list:
  # Make run directory for VPR
  if not os.path.exists(outputDir):
    os.mkdir(outputDir)
  
  # For each file in the directory, run VPR
  for fileName in os.listdir(blifDir):
  
    # Get benchmark name
    fileNameList = fileName.split(".")
    benchmark = str(fileNameList[0])
  
    # Make benchmark directory
    if not os.path.exists(outputDir+"/"+benchmark):
      os.mkdir(outputDir+"/"+benchmark)
  
    # Copy blif from mockturtle to vpr
    copyfile(blifDir+fileName, outputDir+"/"+benchmark+"/"+fileName);
  
    # Run vpr
    value = channelWidthDict[benchmark]
    value = int(value*1.25)
    if value%2 != 0:
      value = value+1
   
    print(benchmark + ":")
    print ("  Trying channel width "+str(value))
    os.system(vprDir+"vpr "+xmlFile+" "+outputDir+"/"+benchmark+"/"+fileName+" "+vprOptions+" "+str(value)+" --seed "+str(seed)+" > /dev/null") 
    with open ('vpr_stdout.log') as vpr_output:
      if 'Routing failed' in vpr_output.read():
        print ("    Routing failed for " + benchmark + " with channel width " + str(value));
        print ("  Trying channel width "+str(value-2))
        os.system(vprDir+"vpr "+xmlFile+" "+outputDir+"/"+benchmark+"/"+fileName+" "+vprOptions+" "+str(value-2)+" --seed "+str(seed)+" > /dev/null") 
    with open ('vpr_stdout.log') as vpr_output:
      if 'Routing failed' in vpr_output.read():
        print ("    Routing failed for " + benchmark + " with channel width " + str(value-2));
        print ("  Trying channel width "+str(value+2))
        os.system(vprDir+"vpr "+xmlFile+" "+outputDir+"/"+benchmark+"/"+fileName+" "+vprOptions+" "+str(value+2)+" --seed "+str(seed)+" > /dev/null") 
    with open ('vpr_stdout.log') as vpr_output:
      if 'Routing failed' in vpr_output.read():
        print ("    Routing failed for " + benchmark + " with channel width " + str(value+2));
        continue;
  
    #for incr in range(0,2): 
    #  os.system(vprDir+"vpr "+xmlFile+" "+outputDir+benchmark+"/"+fileName+" "+vprOptions+" "+str(value+incr*2)) 
    #  with open ('vpr_stdout.log') as vpr_output:
    #    if 'Routing failed' not in vpr_output.read():
    #      break
    #  os.system(vprDir+"vpr "+xmlFile+" "+outputDir+benchmark+"/"+fileName+" "+vprOptions+" "+str(value-(incr+1)*2)) 
    #  with open ('vpr_stdout.log') as vpr_output:
    #    if 'Routing failed' not in vpr_output.read():
    #      break
  
    os.rename(runDir+"vpr_stdout.log", outputDir+"/"+benchmark+"/vpr_stdout.log")
    if os.path.exists(runDir+benchmark+".net"):
      os.rename(runDir+benchmark+".net", outputDir+"/"+benchmark+"/"+benchmark+".net")
    if os.path.exists(runDir+benchmark+".place"):
      os.rename(runDir+benchmark+".place", outputDir+"/"+benchmark+"/"+benchmark+".place")
    if os.path.exists(runDir+benchmark+".route"):
      os.rename(runDir+benchmark+".route", outputDir+"/"+benchmark+"/"+benchmark+".route")
    if os.path.exists(runDir+"report_timing.hold.rpt"):
      os.rename(runDir+"report_timing.hold.rpt", outputDir+"/"+benchmark+"/report_timing.hold.rpt")
    if os.path.exists(runDir+"report_timing.setup.rpt"):
      os.rename(runDir+"report_timing.setup.rpt", outputDir+"/"+benchmark+"/report_timing.setup.rpt")
  
  os.rename(outputDir, outputDir+"_"+str(seed))

os.mkdir(outputDir);
for seed in seed_list:
  os.rename("output_"+str(seed), outputDir+"/output_"+str(seed));
  os.system("./summarize_timing_report.py "+outputDir+"/output_"+str(seed)+" > "+outputDir+"/vpr_result_"+str(seed)+".txt");


