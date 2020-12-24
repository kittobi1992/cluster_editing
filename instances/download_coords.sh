#!/bin/bash

curl "https://cloud.iti.kit.edu/index.php/s/NoWrmMgAWnzKyt3/download?path=/&files=exact_coords.zip" -o "exact_coords.zip"
curl "https://cloud.iti.kit.edu/index.php/s/NoWrmMgAWnzKyt3/download?path=/&files=heur_coords.zip" -o "heur_coords.zip"

unzip exact_coords.zip
unzip heur_coords.zip

rm exact_coords.zip heur_coords.zip 
