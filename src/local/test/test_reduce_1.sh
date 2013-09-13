#!/bin/sh

python ../reduce.py --input reduce_input_1/* --output=reduce_output_1 --name 'TestReduce1' --force --num-tasks 3 --sort-fields 3 --bin-fields 1 --sort-size=10000 -- cat

