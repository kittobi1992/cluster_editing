# TU-Berlin cluEdit Wrapper

Code from: https://fpt.akt.tu-berlin.de/cluEdit/

## How to run stuff

The input format is the same as of the JENA project (`.cm`, see the example
`p4.cm` a path of 4 edges).

A `.cm` file can be passed to the executable as follows:

    java -jar ClusterEditOpt.jar p4.cm NULL 0 0
    
I'm not sure what the last three arguments are doing, but this seems to work.

## How to get `.cm` files

Use the python script `pace2cm.py` to convert PACE format input graphs to the
`.cm` format.

    > ./pace2cm.py --help
    usage: pace2cm.py [-h] in_path cm_path

    positional arguments:
      in_path     PACE format input file
      cm_path     Jena / PEACE format (.cm) output file

    optional arguments:
      -h, --help  show this help message and exit
