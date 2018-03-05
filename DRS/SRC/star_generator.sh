#!/bin/bash

cd /home/main/Documents/DRS/SRC

NBR_ETOILES=1

for ((etoile = 0; etoile < $NBR_ETOILES; etoile++)); do
    python narval_simu2.py st0
done

