#!sh
for i in {1..32} 
do
  R --vanilla --slave -f microbenchmark-local.R --args $i
done
