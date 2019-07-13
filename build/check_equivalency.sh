
FILES=/home/kimjin14/work/mapping/mig_mapping/MIG_experiments/MIG/*

for f in $FILES
do
  file=${f##*/}
  OIFS=$IFS
  IFS='.'
  benchmark_name=(${file})
  IFS=$OIFS
  #echo "Processing $file file..."

  # run benchmark
  #./run_benchmarks/run_benchmark $f

  

  ./run_benchmarks/run_benchmark $f > log/${benchmark_name[1]}.log
 
  diff blif/${benchmark_name[1]}.blif blif/${benchmark_name[1]}_carry.blif >> log/${benchmark_name[1]}.log

  ../../abc/abc -c "cec blif/${benchmark_name[1]}.blif blif/${benchmark_name[1]}_carry.blif" >> log/${benchmark_name[1]}.log

  echo "check log/${benchmark_name[1]}.log for results."

  
  # parse results
  #python parse_log.py output_$file.log 

done

grep "Networks are NOT EQUIVALENT" log/* 
