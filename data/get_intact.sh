#! /usr/bin/bash

## run from the data directory
wget ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip
unzip intact.zip

## As this is large file "intact.txt" of 4.9 GB I will not be submitting that downloaded file

## cut columns we need, only using intact.txt
cut -f 1,2,5,6,7,9,10,11,12,13,14,17,18,19,20 intact.txt > intact_cut.txt


## based on advice from IntAct try to use biological role to determine positive or negative regulation
## originally columns 17, 18 from https://psicquic.github.io/MITAB27Format.html
## NOT ! (both Biological roles = (unspecified role))
## For directionality and/or positive or negative
## This gives us just over 10,000 to work with for direction and whether positive or negative
awk -F '\t' '($12 !~ /(unspecified role)/ || $13 !~ /(unspecified role)/)' intact_cut.txt > intact_cut_direction.txt

## just for associations. Initially going to do 1000, now all future work.
## non-directional in bi-directional querying in neo4j
awk -F '\t' '($12 ~ /(unspecified role)/ && $13 ~ /(unspecified role)/)' intact_cut.txt | head -n1000 > intact_associations.txt
