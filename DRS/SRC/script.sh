#!/bin/bash
################################################################################
# Script pour le cas ou l'ensemble des fichiers est deposé d'un coup dans un   #
# seul dossier.                                                                #
################################################################################


. functions.sh


################################################################################
####### DEPLACEMENT DES FICHIERS DANS DES DOSSIERS PROPRES A CHAQUE NUIT #######
################################################################################

memoire="../Brut_sas/Narval_19950524_171500_fla.fts"
for fich in ../Brut_sas_ala/*ts; do
	
    if [  -e fich ]; then
    echo Ups
	break
    fi
    
    if [[ $fich =~ (Narval|simu[0-9])_[0-9]{8}_[0-9]{6}_[a-z0-9]{3}\.fts ]]; then
	fichier=$fich
	echo $fichier
    elif [[ $fich =~ [0-9]+[a-z]\.fits ]]; then
	fichier=$(python renamator.py $fich)
    fi
    
    if (( $(($(to_sec $fichier) - $(to_sec $memoire))) < 72000 )); then
	mv $fichier "../Brut_sas_ala/$(get_date $memoire)"
    else
	mkdir -p "../Brut_sas_ala/$(get_date $fichier)" ;  mv $fichier $_
	memoire=$fichier
    fi
    rm -f $fich
   
done

################################################################################
######################## TRAITEMENT DES FICHIERS ###############################
################################################################################

for dossier in ../Brut_sas_ala/*; do

    echo $dossier
    if [ -f "$dossier" ]; then
	continue
    fi
    if [[ $dossier =~ [0-9]{8}$ ]]; then
	rm -f /home/main/Documents/DRS/Brut/*
	rm -f /home/main/Documents/DRS/FILES/*
	
	echo DEPOUILLEMENT DU DOSSIER = $dossier
	
	mkdir -p $dossier/stars
	/home/main/Documents/DRS/SRC/arthur2.sh &
	sleep 1s
	pid_arthur=$!


	for fichier in $dossier/*; do

	    # Pour chaque fichier dans chaque nuit, on copie les fichiers propres au
	    # pretraitement dans Brut avec arthur activé. La nuit commence...

	    # Les fichiers stars sont copiés quant à eux dans un dossier à part pour
	    # une utilisation future
	    
	  if [[ $fichier =~ bia.fts$ ]] || verif_flat $dossier $fichier || verif_fap $dossier $fichier; then
		cp $fichier /home/main/Documents/DRS/Brut/
		echo fichier preparatoire copie ! $fichier
		sleep 2s
	    elif [[ $fichier =~ (st0|th0)\.fts$ ]]; then
		cp $fichier $dossier/stars
		echo fichier etoile copie ! $fichier
	    fi
	done
        wait $pid_arthur

    
 
	/home/main/Documents/DRS/SRC/stars2.sh &
	
	sleep 10s
	# On lance maintenant le traitement des étoiles
	
	for etoile in ./$dossier/stars/*; do
		echo Depouillement de $etoile
	    cp $etoile /home/main/Documents/DRS/Brut
	    rm $etoile
	    sleep 10s
	done

	sleep 3s

	# On attend que la dernière étoile de la nuit supprime son .star. Et donc
	# cela indique que la nuit est terminée.
	
	last=$(ls -Art /home/main/Documents/DRS/SRC/.star* | tail -n 1)
	inotifywait -qq -e delete $last
	
    fi
done

	    
