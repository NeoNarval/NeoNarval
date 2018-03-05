#!/bin/bash

cd /home/main/Documents/DRS/SRC


NOMBRE_FLAT=1

rm ../Brut/*
rm ../FILES/*

for ((a=0 ; a < $NOMBRE_FLAT ; a++))
do
  python narval_simu2.py fla
  sleep 5s 
done

python narval_simu2.py fp0
sleep 5s

python bias_generator.py
sleep 5s

# python th_generator.py
# sleep 5s

# python narval_simu2.py st0
# sleep 5s

