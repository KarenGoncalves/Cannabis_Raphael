cut -f 1,2,3,5,10 metadata/metadata.txt |
 awk -F "\t" 'BEGIN {OFS="\t"; print "replicateName", "BioProject", "Layout", "Length_type", "tissue", "Part"}
{
  if (tolower($5) ~ /flo/ || tolower($5) ~ /bud/) {tissue="Flower"}
  else if (tolower($5) ~ /meristem/) {tissue="Meristem"}
  else if (tolower($5) ~ /root/) {tissue="Root"}
  else if (tolower($5) ~ /stem/ || tolower($5) ~ /shoot/) {tissue="Stem"}
  else if (tolower($5) ~ /leaf/ || tolower($5) ~ /petiole/) {tissue="Leaf"}
  else {tissue="Trichome"};
  print $1, $2, $3, $4, tissue, $5
} ' > metadata/metadata_pca.txt
