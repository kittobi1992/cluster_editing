#/bin/bash

curl https://fpt.akt.tu-berlin.de/pace2021/exact.tar.gz -o exact.tar.gz
curl https://fpt.akt.tu-berlin.de/pace2021/heur.tar.gz -o heur.tar.gz
tar -xf exact.tar.gz
tar -xf heur.tar.gz
rm exact.tar.gz heur.tar.gz
