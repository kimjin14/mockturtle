#!/usr/bin/env python

import sys 
import os
from os import path

def parse_delay(line): 
  delays = line.split(')');
  if (len(delays) > 1):
    delays = delays[1].split();
    if (len(delays) <= 0):
      return 0.0
    else:
      return delays[0]
  else:
    return 0.0

def parse_file(file_name):

  # Store routing delays
  routing_value = 0;
  routing_count = 0;
  carry_lut_routing_value = 0;
  carry_lut_routing_count = 0;
  carry_5lut_routing_value = 0;
  carry_5lut_routing_count = 0;
  
  # Store LUT delays
  lut_value = 0;
  lut_count = 0;
 
  carry_cin_value = 0;
  carry_cin_count = 0;
  carry_lut_value = 0;
  carry_lut_count = 0;
  carry_5lut_value = 0;
  carry_5lut_count = 0;
 
  flag_start = 0;
  flag_lut = 0;
  flag_carry = 0;
  flag_cin = 0;
  flag_cout = 0;
  
  hold_value = 0;
  flag_hold = 0;
  with open(file_name) as report1:
    for line in report1:
    
      if "Path 1\n" in line:
        flag_start = 1
      elif "Path 2\n" in line:
        flag_start = 0
        break;
      # Process first path only 
      if (flag_start == 1):
        # This indicates LUT mapping only
        if "(.names" in line:
          flag_lut = 1;
          flag_carry = 0;
        # This indicates LUT + carry - separate to through LUT or cout
        elif "(lut_adder" in line:
          flag_carry = 1;
          flag_lut = 0;
          flag_cin = 0;
          flag_sumout = 0;
          flag_out = 0;
          if "cin" in line:
            flag_cin = 1;
          elif "sumout" in line:
            flag_sumout = 1;
          elif "out" in line:
            flag_out = 1;

        delay = parse_delay(line);
        if (flag_lut == 1):
          if "primitive" in line:
            lut_count = lut_count + 1;
            lut_value = lut_value + float(delay);
          elif "routing" in line:
            routing_count = routing_count + 1;
            routing_value = routing_value + float(delay);
        elif (flag_carry == 1):
          if "primitive" in line:
            if (flag_cin == 1):
              carry_cin_count = carry_cin_count + 1;
              carry_cin_value = carry_cin_value + float(delay);
            else:
              flag_hold = 1;
              hold_value = float(delay);
          elif "routing" in line:
            if (flag_sumout == 1):
              if (flag_hold == 1):
                carry_lut_count = carry_lut_count + 1;
                carry_lut_value = carry_lut_value + hold_value;
                flag_hold = 0;
              carry_lut_routing_count = carry_lut_routing_count + 1;
              carry_lut_routing_value = carry_lut_routing_value + float(delay);
              hold_value = 0;
            elif (flag_out == 1):
              if (flag_hold == 1):
                carry_5lut_count = carry_5lut_count + 1;
                carry_5lut_value = carry_5lut_value + hold_value;
                flag_hold = 0;
              carry_5lut_routing_count = carry_5lut_routing_count + 1;
              carry_5lut_routing_value = carry_5lut_routing_value + float(delay);
              hold_value = 0; 
     
  total_delay = (lut_value + routing_value);
  total_delay += (carry_lut_value + carry_lut_routing_value); 
  total_delay += (carry_5lut_value + carry_5lut_routing_value); 
  total_delay += (carry_cin_value);
  

  #print file_name,
  print str(lut_count)+"("+str(lut_value)+"),",
  print str(routing_count)+"("+str(routing_value)+"),",
  print str(carry_lut_count)+"("+str(carry_lut_value)+"),",
  print str(carry_lut_routing_count)+"("+str(carry_lut_routing_value)+"),",
  print str(carry_5lut_count)+"("+str(carry_5lut_value)+"),",
  print str(carry_5lut_routing_count)+"("+str(carry_5lut_routing_value)+"),",
  print str(carry_cin_count)+"("+str(carry_cin_value)+"),",
  print str(lut_value + carry_lut_value + carry_5lut_value + carry_cin_value)+",",
  print str(routing_value + carry_lut_routing_value + carry_5lut_routing_value)+",",
  print str(total_delay)
  
  return total_delay;  

def main():

  rundir = sys.argv[1]
  allsubdir = os.listdir(rundir)

  geomean = 1.0; 
  n = 0;

  for subdir in allsubdir:
    print subdir + "&",
    vpr_file = rundir + '/' + subdir + '/vpr_stdout.log'
    vprpassflag = 1
    if not path.exists(vpr_file):
      print "Increase channel width manually."
      continue;
    with open(vpr_file) as vprout:
      for line in vprout:
        if 'Circuit is unroutable' in line:
          #print line
          vprpassflag = 0
    if vprpassflag is 1: 
      report = rundir + '/' + subdir + '/report_timing.setup.rpt'
      if path.exists(report):
        geomean *= parse_file (report);
        n = n+1;
      else:
        print "Did not place and route correctly."
    else:
      print 'cannot route.'

  print 'geomean '+str(n)+','+str(geomean**(1/float(n)))+'\n';

if __name__ == "__main__":
    main()
