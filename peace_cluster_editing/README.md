This is the source code for PEACE obtained via https://bio.informatik.uni-jena.de/software/peace/. It has been been
modified in order to successfully compile it and `peace.py` has been added as an utility.

## Building

```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make weighted_cluster_editing
cd ..
```


## Run executable

```bash
mkdir -p results
./build/weighted_cluster_editing --unweighted --mode 2 data/c1.cm results/c1.out
```


## Run PEACE on a PACE instance
```bash
./peace.py --binary ./build/weighted_cluster_editing --instance ../instances/exact/exact031.gr
```