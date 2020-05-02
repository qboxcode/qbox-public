#!/bin/bash
#
# qbox_species_temp.sh
#
# compute the temperature of each atom of a given species in an MD run
#
# use: qbox_species_temp.sh species_name file [file ...]
#
if (($#<2))
 then echo " use: qbox_species_temp.sh species file [file ...]"
      exit 1
fi
species=$1
shift 1
echo "# temperature of atoms of species: $species, file: ${*}"
# extract mass of species
xpath='//species[@name="'$species'"]/mass'
elem=$(basename $xpath)
mass=$(xmllint --xpath "$xpath" ${*} |sed "s/<$elem>//g" |sed "s/<\/$elem>/>/g" |tr '>' '\n' | sed "/^$/d")
# Boltzmann constant in hartree/K
kB=3.1667907e-06
# nucleon mass
nmass=1822.89
# conversion factor from velocity^2 (a.u.) to K
# fac = (2/3) * (1/2) * nmass * mass * (1/kB)
# extract velocity of each atom of the given species
xpath='//atomset/atom[@species="'$species'"]/velocity'
elem=$(basename $xpath)
xmllint --xpath "$xpath" ${*} |sed "s/<$elem>//g" |sed "s/<\/$elem>/>/g" |tr '>' '\n' | sed "/^$/d" | \
awk -v kB=$kB -v nmass=$nmass -v mass=$mass \
'{v2=$1*$1+$2*$2+$3*$3; print (2.0/3.0)*0.5*nmass*mass*v2/kB}' -
