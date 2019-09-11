
FILES=/home/kimjin14/work/mapping/mig_mapping/MIG_experiments/MIG/*

for f in $FILES
do
  file=${f##*/}
  OIFS=$IFS
  IFS='.'
  benchmark_name=(${file})
  IFS=$OIFS
  echo "Processing $file file..."

  # run mapping on benchmark
  ./run_benchmarks/run_benchmark $f > log/${benchmark_name[1]}.log
  diff blif/${benchmark_name[1]}.blif blif/${benchmark_name[1]}_carry.blif >> log/${benchmark_name[1]}.log
  ../../abc/abc -c "cec blif/${benchmark_name[1]}.blif blif/${benchmark_name[1]}_carry.blif" >> log/${benchmark_name[1]}.log

  # check if mapping produced correct blif
  echo "check log/${benchmark_name[1]}.log for results."
  grep "Networks are equivalent" log/${benchmark_name[1]}.log

done

grep "Networks are NOT EQUIVALENT" log/* 
