#!/bin/bash
# density = 0.01
python generate_random.py --m regular --n 200 --d 2 --seed 42 --output /home/zgh23/code/SparseWS/data/origin
python generate_random.py --m regular --n 600 --d 6 --seed 42 --output /home/zgh23/code/SparseWS/data/origin
python generate_random.py --m regular --n 1000 --d 10 --seed 42 --output /home/zgh23/code/SparseWS/data/origin
python generate_random.py --m regular --n 3000 --d 30 --seed 42 --output /home/zgh23/code/SparseWS/data/origin
python generate_random.py --m regular --n 5000 --d 50 --seed 42 --output /home/zgh23/code/SparseWS/data/origin
python generate_random.py --m regular --n 7000 --d 70 --seed 42 --output /home/zgh23/code/SparseWS/data/origin
python generate_random.py --m regular --n 10000 --d 100 --seed 42 --output /home/zgh23/code/SparseWS/data/origin
