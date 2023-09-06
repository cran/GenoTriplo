#!/bin/bash
# PLEASE CLOSE THIS FILE
fichier=$1
oldIFS=$IFS
IFS=$'\n'
# PLEASE CLOSE THIS FILE

# If This open, close it ! It is because you don't have the possibility to launch bash script.
# Don't worry, the creation of the dataset will proceed once you close this window !

# _____  _      ______           _____ ______      _____ _      ____   _____ ______ 
#|  __ \| |    |  ____|   /\    / ____|  ____|    / ____| |    / __ \ / ____|  ____|
#| |__) | |    | |__     /  \  | (___ | |__      | |    | |   | |  | | (___ | |__   
#|  ___/| |    |  __|   / /\ \  \___ \|  __|     | |    | |   | |  | |\___ \|  __|  
#| |    | |____| |____ / ____ \ ____) | |____    | |____| |___| |__| |____) | |____ 
#|_|    |______|______/_/    \_\_____/|______|    \_____|______\____/|_____/|______|
#                                                                                 
#  _______  _    _  _____   _____     ______  _____  _       ______ 
# |__   __|| |  | ||_   _| / ____|   |  ____||_   _|| |     |  ____|
#    | |   | |__| |  | |  | (___     | |__     | |  | |     | |__   
#    | |   |  __  |  | |   \___ \    |  __|    | |  | |     |  __|  
#    | |   | |  | | _| |_  ____) |   | |      _| |_ | |____ | |____ 
#    |_|   |_|  |_||_____||_____/    |_|     |_____||______||______|
#










pattern1="_[A-Z][0-9][0-9].CEL"
pattern2="_[A-Z][0-9].CEL"

nbL=$(cat $fichier | wc -l)

var=$(grep -n -w '^probeset_id' $fichier | cut -d ":" -f 1)

# Selectionne les lignes a partir de 'probeset_id'
# Change le pattern1 par rien (il y a rien en les deux //) puis le pattern2
(sed -n $var,$nbL'p' $fichier | sed "s/$pattern1//g" | sed "s/$pattern2//g") > './tmp.txt'

