from shutil import copyfile
import os
import sys

channelWidthDict = {
  "aes_core":72, "comp":190, "CSA464":62, "des_area":88, "diffeq1":98, "div16":98, "hamming":74, "i2c":44, "log2":170, "MAC32":100, "max":74, "mem_ctrl":104, "MUL32":104, "mult64":132, "pci_spoci_ctrl":64, "rca_128":158, "rca_256":66, "rca_32":62, "rca_64":58, "revx":78, "sasc":48, "simple_spi":40, "spi":76, "sqrt32":78, "square":118, "ss_pcm":40, "systemcaes":76, "systemcdes":50, "tern_add":66, "tv80":110, "usb_phy":38,
  "aes_core_carry":72, "comp_carry":190, "CSA464_carry":62, "des_area_carry":88, "diffeq1_carry":98, "div16_carry":98, "hamming_carry":74, "i2c_carry":44, "log2_carry":170, "MAC32_carry":100, "max_carry":74, "mem_ctrl_carry":104, "MUL32_carry":104, "mult64_carry":132, "pci_spoci_ctrl_carry":64, "rca_128_carry":158, "rca_256_carry":66, "rca_32_carry":62, "rca_64_carry":58, "revx_carry":78, "sasc_carry":48, "simple_spi_carry":40, "spi_carry":76, "sqrt32_carry":78, "square_carry":118, "ss_pcm_carry":40, "systemcaes_carry":76, "systemcdes_carry":50, "tern_add_carry":66, "tv80_carry":110, "usb_phy_carry":38
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

