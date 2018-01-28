#!sh
for i in {13..48}
do
  R --vanilla --slave -f microbenchmark-local.R --args $i
done
