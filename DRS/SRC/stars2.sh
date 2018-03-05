#!/bin/bash

################################################################################
# Script gérant les arrivées d'étoiles et de thorium-argon au cours de la nuit #
################################################################################

. functions.sh
NB_STAR=1			# Nombre de fichiers d'étoiles nécessaires pour une étoile
STAR_DONE=0

NB_TH=1				# Nombre de thorium-argon nécessaires pour un thorium-argon
TH_DONE=0

################################################################################
######################## VERIFICATION DES LOCK #################################
################################################################################

if [ ! -e  ".reussi" ]; then
    echo "Fichier reussi absent, le script s'arrête"
    exit 1
fi

if [ -e ".lockstar" ]
then
    echo "Fichier lock present, le script s'arrête"
    exit 1
fi

echo "lock file star" >> .lockstar

trap "rm .lockstar" EXIT	# Au cas ou le programme plante ou est killed, on supprime le lock

echo "\
################################################################################
#################### ATTENTE DES ETOILES DURANT LA NUIT ########################
################################################################################"

while [ true ]; do
    FICHIER="../Brut/"$(inotifywait -q -e create -e moved_to --format %f ../Brut/.)

    if [[ $FICHIER =~ st0\.fts$ ]]; then
	echo "stars2: Star found and queued!"
	stars[$STAR_DONE]=$FICHIER
	(( STAR_DONE++ ))
    elif [[ $FICHIER =~ th0\.fts$ ]]; then
	echo "stars2: TH found and queued!"
	th[$TH_DONE]=$FICHIER
	(( TH_DONE++ ))
    fi

    if [[ $STAR_DONE == $NB_STAR ]]; then
	STAR_DONE=0
	echo "stars2: Script stars3:stars lancé !"
	./stars3.sh st ${stars[@]} & 
    elif [[ $TH_DONE == $NB_TH ]]; then
	TH_DONE=0
	echo "stars2: Script stars3:thorium lancé !"
	./stars3.sh th  ${th[@]} &
    fi   
    
done

