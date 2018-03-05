#!/bin/bash

cd /home/main/Documents/DRS/SRC/
VOIES=2

set_param ()
{
    sed -i.bak "s~^\($1 : \).*~\1$2~" $3
}

get_param ()
{
    awk -F" : " "/^$1/"' {print $2}' $2
}

add_to_list ()
{
    sed -i.bak "/^\[$1\]/ { N; s~~\[$1\]\n$2~ }" $3
}

A_name ()
{
    echo Amatrix$1\_OR$2\_LA$3.npz
}

initialize ()
{
    echo -e "[BIAS]\n\n[FLAT]\n\n[FABRY]" > $1
}

to_sec ()
{
    local fichier=$1
    local date=$(echo $fichier | grep -oE "[0-9]{8}")
    local heure=$(echo $fichier | grep -oE "_[0-9]{6}_")
    local fulldate="$date ${heure:1:2}:${heure:3:2}:${heure:5:2}"
    local sec=$(date --universal -d "$fulldate" +%s)
    echo $sec
}

get_date ()
{
    local fichier=$1
    local date=$(echo $fichier | grep -oE "_[0-9]{8}_")
    echo ${date:1:8}
}


# Permet de verifier qu'un flat est bien du début de la nuit.
# Inutile en véritable utilisation mais évite de relancer envelope.py trop de fois
# lors de l'utilisation en mode script.sh.

verif_flat ()
{
    local dossier=$(echo $1 | grep -oE "[0-9]{8}")
    local fichier=$2
    if [[ $fichier =~ fla\.fts$ ]]; then
	if (( $(get_date $fichier) == $dossier )); then
	    return 0
	else
	    return 1
	fi
    else
	return 1
    fi
    
}

# Idem pour le Fabry-Pérot ! 
verif_fap ()
{
    local dossier=$(echo $1 | grep -oE "[0-9]{8}")
    local fichier=$2
    if [[ $fichier =~ fp0\.fts$ ]]; then
	if (( $(get_date $fichier) == $dossier )); then
	    return 0
	else
	    return 1
	fi
    else
	return 1
    fi
    
}

get_spectrum ()
{

    local DATA_FILE="../DATA/Amatrix_DATA.txt"
    
    local NBR_ORDRE_TOTAL=40
    local SPECTRE_DONE=0
    local fts_file=$1

    set_param "Test fts file" $fts_file $DATA_FILE
    
    for ((ordre = 0; ordre < $NBR_ORDRE_TOTAL; ordre++)); do
	set_param "Order n°" $ordre $DATA_FILE
	for ((voie = 1; voie <= $VOIES; voie++)); do
	    lvl=$(get_param "Evolution level" $DATA_FILE)
	    set_param "Lane n°" $voie $DATA_FILE
	    set_param "A matrix file" "../TEMP/Amatrix_v$lvl\_OR$ordre\_LA$voie.npz" $DATA_FILE
	    python generate_lane_spectrum.py
	    if [[ $? == 0 ]]; then
		(( SPECTRE_DONE++ ))
        otru=$((ordre+20))
		echo "get_spectrum: Spectre ordre n°$otru, voie n°$voie fait"
	    fi
	done
    done

}

create_A_matrix ()
{

    local NBR_ORDRE_TOTAL=40
    local MATRIX_DONE=0
    local fp_file=$1
    
    for ((ordre = 0; ordre < $NBR_ORDRE_TOTAL; ordre++)); do
	set_param "Order n°" $ordre ../DATA/Amatrix_DATA.txt
	for ((voie = 1; voie <= $VOIES; voie++)); do
	    set_param "Lane n°" $voie ../DATA/Amatrix_DATA.txt
            python create_A_matrix.py
	    if [[ $? == 0 ]]; then    
		((MATRIX_DONE++))
	    fi
            echo $MATRIX_DONE "matrices de 80 faites par arthur2 et create_A_matrix"
	done
    done
}

cosmic ()
{
    # ARGUMENT 1 : fichier à traiter
    # ARUGMENT 2 : fichier de sortie
    
    local NBR_ORDRE_TOTAL=40
    local fichier_fts=$(basename $1)
    local lvl=$(get_param "Evolution level" ../DATA/Amatrix_DATA.txt)
    echo $fichier_fts
    
    for ((ordre = 0; ordre < $NBR_ORDRE_TOTAL; ordre++)); do
	set_param "Order n°" $ordre ../DATA/Amatrix_DATA.txt
	for ((voie = 1; voie <= $VOIES; voie++)); do
	    set_param "Lane n°" $voie ../DATA/Amatrix_DATA.txt
	    set_param "A matrix file" "../TEMP/Amatrix_v${lvl}_OR${ordre}_LA$voie.npz" ../DATA/Amatrix_DATA.txt
	    python cosmic_auto.py $2
	done
    done
}


