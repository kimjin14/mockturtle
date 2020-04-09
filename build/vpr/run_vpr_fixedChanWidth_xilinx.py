#!/usr/bin/env python3

from shutil import copyfile
import os
import sys

channelWidthDict = {
  "aes_core":62, "comp":178, "CSA464":62, "des_area":80, "diffeq1":76, "div16":86, "hamming":76, "i2c":42, "log2":162, "MAC32":90, "max":62, "mem_ctrl":88, "MUL32":94, "mult64":124, "pci_spoci_ctrl":70, "rca_128":148, "rca_256":70, "rca_32":62, "rca_64":56, "revx":70, "sasc":38, "simple_spi":40, "spi":66, "sqrt32":76, "square":110, "ss_pcm":36, "systemcaes":86, "systemcdes":50, "tern_add":70, "tv80":98, "usb_phy":38,
  "aes_core_carry":62, "comp_carry":178, "CSA464_carry":62, "des_area_carry":80, "diffeq1_carry":76, "div16_carry":86, "hamming_carry":76, "i2c_carry":42, "log2_carry":162, "MAC32_carry":90, "max_carry":62, "mem_ctrl_carry":88, "MUL32_carry":94, "mult64_carry":124, "pci_spoci_ctrl_carry":70, "rca_128_carry":148, "rca_256_carry":70, "rca_32_carry":62, "rca_64_carry":56, "revx_carry":70, "sasc_carry":38, "simple_spi_carry":40, "spi_carry":66, "sqrt32_carry":76, "square_carry":110, "ss_pcm_carry":36, "systemcaes_carry":86, "systemcdes_carry":50, "tern_add_carry":70, "tv80_carry":98, "usb_phy_carry":38
}

blifDir= "/home/kimjin14/work/mapping/mig_mapping/mockturtle/build/"
vprDir="/home/kimjin14/work/vtr/vtr-verilog-to-routing/vpr/"
runDir="/home/kimjin14/work/mapping/mig_mapping/mockturtle/build/vpr/"
outputDir="/home/kimjin14/work/mapping/mig_mapping/mockturtle/build/vpr/output/"

xmlFile = "arch/Xilinx_style_arch.xml"
vprOptions = "--timing_report_detail aggregated --route_chan_width"

# Get directory of what you want to run
blifDir=blifDir+sys.argv[1]+"/"

# Make run directory for VPR
if not os.path.exists(outputDir):
  os.mkdir(outputDir)

# For each file in the directory, run VPR
for fileName in os.listdir(blifDir):

  # Get benchmark name
  fileNameList = fileName.split(".")
  benchmark = str(fileNameList[0])

  # Make benchmark directory
  if not os.path.exists(outputDir+benchmark):
    os.mkdir(outputDir+benchmark)

  # Copy blif from mockturtle to vpr
  copyfile(blifDir+fileName, outputDir+benchmark+"/"+fileName);

  # Run vpr
  value = channelWidthDict[benchmark]
  value = int(value*1.25)
  if value%2 != 0:
    value = value+1
  os.system(vprDir+"vpr "+xmlFile+" "+outputDir+benchmark+"/"+fileName+" "+vprOptions+" "+str(value)) 

  print(benchmark)
  os.rename(runDir+"vpr_stdout.log", outputDir+benchmark+"/vpr_stdout.log")
  if os.path.exists(runDir+benchmark+".net"):
    os.rename(runDir+benchmark+".net", outputDir+benchmark+"/"+benchmark+".net")
  if os.path.exists(runDir+benchmark+".place"):
    os.rename(runDir+benchmark+".place", outputDir+benchmark+"/"+benchmark+".place")
  if os.path.exists(runDir+benchmark+".route"):
    os.rename(runDir+benchmark+".route", outputDir+benchmark+"/"+benchmark+".route")
  if os.path.exists(runDir+"report_timing.hold.rpt"):
    os.rename(runDir+"report_timing.hold.rpt", outputDir+benchmark+"/report_timing.hold.rpt")
  if os.path.exists(runDir+"report_timing.setup.rpt"):
    os.rename(runDir+"report_timing.setup.rpt", outputDir+benchmark+"/report_timing.setup.rpt")

