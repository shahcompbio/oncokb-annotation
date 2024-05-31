for patient in 107C1 136PC2; do 
    echo "[`date`] Run ANNOVAR for $patient"
    for avinput in `\ls ../results/avinput/${patient}.*.avinput`; do 
        annotate_variation.pl -geneanno -buildver hg38 $avinput ~/data1/packages/annovar/humandb/
    done
done

