BEGIN	{FS=",";}
NR==1 {print;} 
NR>1 && $1 != "NA" {a[i++]=$0;} 
NR > 1 && $1 == "NA" {if (i>=expected_motif_length) {if (c>0) print "NA,NA"; for(b=0;b<i;b++) print a[b]; c++;} i=0;} 
END {if (i>=expected_motif_length) {if (c>0) print "NA,NA"; for(b=0;b<i;b++) print a[b];}}