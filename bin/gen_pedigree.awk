#!/usr/bin/awk -f
#pedigree file structure:
 #       Family ID
 #       Individual ID
 #       Paternal ID
 #       Maternal ID
 #       Sex (1=male; 2=female; other=unknown)
 #       Phenotype#
#
 #       Affected status should be coded as follows:
 #       -9 missing
 #       0 missing
 #       1 unaffected
 #       2 affected

BEGIN {
    OFS="\t"
}  # Begin section

{
    gender=0;
    if ($2 == "xy" || $2 == "XY")
	gender = 1
    else if ($2 == "xx" || $2 == "XX")
	gender = 2
    phenotype = $3 + 1
    
    print $1, $4, 0, 0, gender, phenotype
}        # Loop section
END{}     # End section
