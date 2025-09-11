GENOME="Mada106.fna"

gapseq find -p all "$GENOME"

#FIND TRANSPORTER
gapseq find-transport "$GENOME"

#CREATING DRAFT MODEL
gapseq draft -r "Mada106-all-Reactions.tbl" \
             -t "Mada106-Transporter.tbl" \
             -p "Mada106-all-Pathways.tbl" \
             -c "$GENOME"

#GAP FILLING
#NO SV FILE ATTACHED AS I DON'T KNOW WHAT IT PROVIDES BUT IT SHOULD STILL WORK
gapseq fill -m "Mada106-draft.RDS" \
            -c "Mada106-rxnWeights.RDS" \
            -g "Mada106-rxnXgenes.RDS"
