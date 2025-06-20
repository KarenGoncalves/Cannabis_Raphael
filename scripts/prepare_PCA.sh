#!/bin/sh

# Select columns from metadata to use in R for coexpression analysis and PCA

cut -f 1,2,3,5,10 metadata/metadata.txt |
 awk -F "\t" 'BEGIN {OFS="\t"; print "replicateName", "BioProject", "Layout", "Length_type", "tissue", "Part"}
{
  if (tolower($5) ~ /flo/ || tolower($5) ~ /bud/) {tissue="Flower"} # any tissue name containing floral, flower or bud  was considered flower
  else if (tolower($5) ~ /meristem/) {tissue="Meristem"}
  else if (tolower($5) ~ /root/) {tissue="Root"}
  else if (tolower($5) ~ /stem/ || tolower($5) ~ /shoot/) {tissue="Stem"} # both stem and shoot were considered stem 
  else if (tolower($5) ~ /leaf/ || tolower($5) ~ /petiole/) {tissue="Leaf"} # petiole samples were considered leaf
  else {tissue="Trichome"};
  print $1, $2, $3, $4, tissue, $5
} ' > metadata/metadata_pca.txt
