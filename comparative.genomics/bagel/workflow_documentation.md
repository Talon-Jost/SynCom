What I did for this workflow was as follows:



\#--------------------------------------------------

\# Description and Background for the use of BAGEL5

\#--------------------------------------------------



I ran both the Mada.85 and Mada.106 genomes through BAGEL5 on their online portal: http://bagel5.molgenrug.nl/



DESCRIPTION: BAGEL4 is a web server that enables users to identify and visualize gene clusters in prokaryotic DNA involved in the biosynthesis of Ribosomally synthesized and Post translationally modified Peptides (RiPPs) and (unmodified) bacteriocins.



\#the new BAGEL5 is not yet published but what is available on their website



Citation:

\[BAGEL4: a user-friendly web server to thoroughly mine RiPPs and bacteriocins

Auke J van Heel\*, Anne de Jong\*, Chunxu Song, Jakob H Viel, Jan Kok, Oscar P Kuipers

\*joint first author

Nucleic Acids Research, Volume 46, Issue W1, 2 July 2018, Pages W278â€“W281

https://doi.org/10.1093/nar/gky383]



this resulted in several outputs which are located in the respective files of 85Bagel and 106Bagel. 

1. An FAA file - Containing the protein and amino acid sequences for each of their classifications.
2. An FNA file - The nucleotide sequence provided from the original FASTA input. 
3. A GBK file - A GenBank file which provides more detail on the classification on the actual classifications themselves including the descriptions
4. A Genes file - containing the different sequence length and their assignment the orf number. 
5. A genetable file - containing a more tabular format which combines this information listed above.
6. A transterm file - Transterm facilitates studies of messenger RNAs and translational control signals. Each messenger RNA (mRNA) from GenBank is extracted and broken into its functional components, its coding sequence, initiation context, termination context, flanking sequence representing its 5′ UTR (untranslated region), 3′ UTR and translational signals. In addition, numerical parameters characterising each coding region in Transterm, including codon and GC bias, are available.





\#---------------------------

\# Pulling all of this into R

\#---------------------------



The R workflow has its own documentation in the actual file itself named: bagel\_comparison (located in this folder with the bagel outputs)



However, briefly - my intention was to pull the information provided by bagel into a singular dataframe to work with and to also compare the outputs of the RiPPs and the bacteriocins that the strains may be producing. This resulted in the csv file named: bagel\_combined\_results.csv 



This workflow was very simple and short - I just pulled the information into the script, combined the tables, and calculated the percent similarity between the protein sequences and between the DNA sequences. I dropped anything that was a 100% match, as it doesn't seem useful to compare proteins which have identical structure as my assumption is they will have the exact same use - which doesn't seem too helpful.





\#---------------------------

\# Protein BLAST \& Nucleotide BLAST

\#---------------------------



Running Protein BLAST



Using the resulting table, it had several annotations which were either blank or had no known categorization. Therefore, what I did was I took the amino acid sequence assigned to those proteins and used BLAST to see if it would come up with anything inherently useful. For several of them, it did. As the results came back with some Janibacter proteins. Others weren't so successful. For the ones that didn't have any success with BLAST, I moved on to the next step with them.



Running Nucleotide BLAST



Using the Nucleotide sequence found in the bagel\_combined\_results.csv file, I took one of the amino acid sequences and ran it through BLAST. What I found wasn't too informative. It came back to Janibacter species with the following alignments:



Description			                               | Scientific Name        |Max Sc. |Tot Sc. |Cover  |E value      |Per.ident   |Acc. Len  |Accession  

Janibacter sp. YB324 chromosome, complete genome	       | Janibacter sp. YB324	|183	 |183	  |100%	  |3.00E-42     |100	     |3369845   |0

Janibacter terrae strain COS04-44 chromosome, complete genome  | Janibacter terrae	|183	 |183	  |100%	  |3.00E-42     |100	     |3394906   |0

Janibacter sp. G349 chromosome, complete genome		       | Janibacter sp. G349	|183	 |183	  |100%	  |3.00E-42     |100	     |3507973	|0

Janibacter sp. CX7 chromosome, complete genome		       | Janibacter sp. CX7	|183	 |183	  |100%	  |3.00E-42     |100	     |3229343	|0

Janibacter terrae strain Sec5.9 chromosome, complete genome    | Janibacter terrae	|183	 |183	  |100%	  |3.00E-42     |100	     |3492674	|0



Truly, this is relatively helpful and unhelpful as it does tell us that this is most certainly belonging to a Janibacter species (yay!) and that means it's relatively specific to the Janibacter strains. The bad news is, it doesn't tell us any real information about the protein itself. Which means there's room for some classification about the protein itself if we wanted to go down this road to try and classify certain proteins which might be inhibitory. This also gives us some ability to do some knockouts that might be useful in determining that, also if that is the road we want to go down.



In routes we might be able to take and people who might be of some help. Joshua Kellog and Carly Schissel do some work with molecule conformation as well as some synthesis of these potential molecules. If we can get these synthesized then we might be able to do some specific testing. Although it might be useful to do the meta-transcriptomics first to see when these proteins might actually be expressed before we run down the synthesis route.





\#----------------------------------------------

\# Interpro: Classification of Protein Families

\#----------------------------------------------



I am not entirely sure what the name of this program is; however, this is what it says on their website at least: \[found here: https://www.ebi.ac.uk/interpro/]



CITATION: Blum M, Andreeva A, Florentino LC, Chuguransky SR, Grego T, Hobbs E, Pinto BL, Orr A, Paysan-Lafosse T, Ponamareva I, Salazar GA, Bordin N, Bork P, Bridge A, Colwell L, Gough J, Haft DH, Letunic I, Llinares-López F, Marchler-Bauer A, Meng-Papaxanthos L, Mi H, Natale DA, Orengo CA, Pandurangan AP, Piovesan D, Rivoire C, Sigrist CJA, Thanki N, Thibaud-Nissen F, Thomas PD, Tosatto SCE, Wu CH, Bateman A.

InterPro: the protein sequence classification resource in 2025

Nucleic Acids Research. 2024, doi: 10.1093/nar/gkae1082



DESCRIPTION: InterPro provides functional analysis of proteins by classifying them into families and predicting domains and important sites. To classify proteins in this way, InterPro uses predictive models, known as signatures, provided by several different databases (referred to as member databases) that make up the InterPro consortium. We combine protein signatures from these member databases into a single searchable resource, capitalising on their individual strengths to produce a powerful integrated database and diagnostic tool.



This was easy enough as it just used the bagel\_combined\_results.csv and used the amino acid sequence to run the classification for each of the remaining one or two sequences which didn't match any function to either BAGEL5 or to BLASTp. I did this with two of the sequences - one that had a classification, and the one that didn't (just to test it out). I found a good amount of information returned by this program which include the following:



Cpfam\_protein6 is the sequence which didn't have a known classification based on the results from BAGEL5, BLASTp, and BLASTn.



What was produced from this program:



1. The DNA sequence used.
2. a feature ontology with much more information in a .gff format
3. a JSON file
4. another Microsoft edge file
5. a file that has something in it that I don't quite know its utility or classification. just that it's a TSV file.















































