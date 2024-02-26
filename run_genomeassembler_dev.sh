#!/bin/bash

# preprocess files
cd ./01_Preprocess/scripts
bash submit.sh
cd ../../

# run proof of concept studies
cd ./scripts/
bash submit.sh
cd ../