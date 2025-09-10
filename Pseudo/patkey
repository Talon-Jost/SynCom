GENOME="Andres.fna"

gapseq find -p all "$GENOME"

#FIND TRANSPORTER
gapseq find-transport "$GENOME"

#CREATING DRAFT MODEL
gapseq draft -r "Andres-all-Reactions.tbl" \
             -t "Andres-Transporter.tbl" \
             -p "Andres-all-Pathways.tbl" \
             -c "$GENOME"

#GAP FILLING
#NO SV FILE ATTACHED AS I DON'T KNOW WHAT IT PROVIDES BUT IT SHOULD STILL WORK
gapseq fill -m "Andres-draft.RDS" \
            -c "Andres-rxnWeights.RDS" \
            -g "Andres-rxnXgenes.RDS"
