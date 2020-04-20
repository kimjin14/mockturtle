from shutil import copyfile
import os
import sys

channelWidthDict = {
  "aes_core":70,
  "comp":192,
  "CSA464_mig":50,
  "CSA464":60,
  "des_area":86,
  "diffeq1":90,
  "div16":102,
  "hamming":74,
  "i2c":46,
  "log2":172,
  "MAC32":100,
  "max":68,
  "mem_ctrl":88,
  "MUL32":104,
  "mult64":132,
  "pci_spoci_ctrl":68,
  "rca_128_mig":48,
  "rca_128":156,
  "rca_256_mig":50,
  "rca_256":74,
  "rca_32_mig":44,
  "rca_32":64,
  "rca_64_mig":46,
  "rca_64":60,
  "revx":78,
  "sasc":42,
  "simple_spi":42,
  "spi":88,
  "sqrt32":78,
  "square":116,
  "ss_pcm":38,
  "systemcaes":86,
  "systemcdes":44,
  "tern_add_mig":62,
  "tern_add":70,
  "tv80":110,
  "usb_phy":40,
  "aes_core_carry":70,
  "comp_carry":192,
  "CSA464_mig_carry":50,
  "CSA464_carry":60,
  "des_area_carry":86,
  "diffeq1_carry":90,
  "div16_carry":102,
  "hamming_carry":74,
  "i2c_carry":46,
  "log2_carry":172,
  "MAC32_carry":100,
  "max_carry":68,
  "mem_ctrl_carry":88,
  "MUL32_carry":104,
  "mult64_carry":132,
  "pci_spoci_ctrl_carry":68,
  "rca_128_mig_carry":48,
  "rca_128_carry":156,
  "rca_256_mig_carry":50,
  "rca_256_carry":74,
  "rca_32_mig_carry":44,
  "rca_32_carry":64,
  "rca_64_mig_carry":46,
  "rca_64_carry":60,
  "revx_carry":78,
  "sasc_carry":42,
  "simple_spi_carry":42,
  "spi_carry":88,
  "sqrt32_carry":78,
  "square_carry":116,
  "ss_pcm_carry":38,
  "systemcaes_carry":86,
  "systemcdes_carry":44,
  "tern_add_mig_carry":62,
  "tern_add_carry":70,
  "tv80_carry":110,
  "usb_phy_carry":40
}

blifDir= "/home/kimjin14/work/mapping/mig_mapping/mockturtle/build/"
vprDir="/home/kimjin14/work/vtr/vtr-verilog-to-routing/vpr/"
runDir="/home/kimjin14/work/mapping/mig_mapping/mockturtle/build/vpr/"

xmlFile = "arch/Intel_style_arch.xml"
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

#  "aes_core":162,        "aes_core_carry":162,          
#  "comp":374,            "comp_carry":374,
#  "des_area":214,        "des_area_carry":214,
#  "diffeq1":218,         "diffeq1_carry":218,
#  "div16":258,           "div16_carry":258,
#  "hamming":190,         "hamming_carry":190,
#  "i2c":82,              "i2c_carry":82,
#  "log2":322,            "log2_carry":322,
#  "MAC32":224,           "MAC32_carry":224,
#  "max":190,             "max_carry":190,
#  "mem_ctrl":216,        "mem_ctrl_carry":216,
#  "MUL32":260,           "MUL32_carry":260,
#  "mult64":304,          "mult64_carry":304,
#  "pci_spoci_ctrl":132,  "pci_spoci_ctrl_carry":132,
#  "revx":198,            "revx_carry":198,
#  "sasc":56,             "sasc_carry":56,
#  "simple_spi":94,       "simple_spi_carry":94,
#  "spi":216,             "spi_carry":216,
#  "sqrt32":176,          "sqrt32_carry":176,
#  "square":268,          "square_carry":268,
#  "ss_pcm":48,           "ss_pcm_carry":48,
#  "systemcaes":204,      "systemcaes_carry":204,
#  "systemcdes":162,      "systemcdes_carry":162,
#  "tv80":268,            "tv80_carry":268,
#  "usb_phy":50,          "usb_phy_carry":50,         
#
#  "CSA464_mig":50,       "CSA464_mig_carry":50,          
#  "CSA464":60,           "CSA464_carry":60,           
#  "rca_128_mig":48,      "rca_128_mig_carry":48,      
#  "rca_128":156,         "rca_128_carry":156,         
#  "rca_256_mig":50,      "rca_256_mig_carry":50,      
#  "rca_256":74,          "rca_256_carry":74,          
#  "rca_32_mig":44,       "rca_32_mig_carry":44,       
#  "rca_32":64,           "rca_32_carry":64,           
#  "rca_64_mig":46,       "rca_64_mig_carry":46,       
#  "rca_64":60,           "rca_64_carry":60,           
#  "tern_add_mig":62,     "tern_add_mig_carry":62,     
#  "tern_add":70,         "tern_add_carry":70          

