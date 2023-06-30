#!/bin/bash
fichier=$1
oldIFS=$IFS
IFS=$'\n'

pattern1="_[A-Z][0-9][0-9].CEL"
pattern2="_[A-Z][0-9].CEL"

nbL=$(cat $fichier | wc -l)

var=$(grep -n -w '^probeset_id' $fichier | cut -d ":" -f 1)

# Selectionne les lignes a partir de 'probeset_id'
# Change le pattern1 par rien (il y a rien en les deux //) puis le pattern2
(sed -n $var,$nbL'p' $fichier | sed "s/$pattern1//g" | sed "s/$pattern2//g") > './tmp.txt'

