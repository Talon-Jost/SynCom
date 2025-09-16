
#vsearch --cluster_fast /Users/mcbletz/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Bletz_Lab/Grad_Students/Talon/MadagascarDataProcessing/AmphiBacMada/AmphibBacMADA_FullDatabase_2023.2r.fasta --uc /Users/mcbletz/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Bletz_Lab/Grad_Students/Talon/MadagascarDataProcessing/AmphiBacMada/AmphibBacMADA_FullDatabase_2023.2r_clusters100.txt --id 1 --strand both

#(vsearch) mcbletz@Mollys-MacBook-Pro-2 CultureData % vsearch --cluster_fast /Users/mcbletz/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Bletz_Lab/Computational_Resources/Datasets/CultureData/MadaCultureDB_seqs.fasta --uc /Users/mcbletz/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/Bletz_Lab/Computational_Resources/Datasets/CultureData/MadaCultureDB_seqs_clustered_per100.txt --id 1 --strand both
##*** MODIFIED HEADERS AND COLUMNS MANUALLY IN EXCEL

## merging Madagascar database metadata (culture info, etc) with Vsearch 100% clusterIDs



## sorted by seq plate, toxonomy string then isolate ID

MadaDBdata = read_tsv("./Computational_Resources/Datasets/CultureData/MadaCultureDB_Metadata.txt")

## This deals with the duplicate entries for the isolates; its not perfect, it just selects the first entry from the sorted file (see above)

MadaDBdata_filt = MadaDBdata %>%
  group_by(Isolate.ID) %>% slice_head() %>% # filters to 1 sequence per group
  ungroup()

## for fasta creation
MadaDBdata_filt_fasta = MadaDBdata_filt %>%
  select(Isolate.ID,Sequence) %>% # selecting the two columns needed for fasta creation
  drop_na(Sequence) # get rid of isolates we no sequence

colnames(MadaDBdata_filt_fasta) <- c("seq.name","seq.text") # headers needed for phylotools 

## Make fasta file
phylotools::dat2fasta(MadaDBdata_filt_fasta,"./Computational_Resources/Datasets/CultureData/MadaCultureDB_seqs_dupsfilt.fasta")

## ex
write_tsv(MadaDBdata_filt_fasta,"./MadaCultureDB_dupsfilt_prefile.txt")


## Vsearch done with fasta from MadaDB file (includes dublicates); this has more isolates than AmphiBac; This is the vsearch output (code above)

MadaDB_VsearchClusters = read_tsv("./Computational_Resources/Datasets/CultureData/MadaCultureDB_seqs_clustered_per100.txt")

## merging clustering with Mada DB metadata that has had the dubs filtered
## NOTE Vsearch was not done with the duplicate filtered fasta...I think that's OK for now.but it does mean that there could be isolates with same ID that are "different" clusters because there are some that didn't ID to the same taxonomy.
MadaDB_VsearchClusters_wdata = left_join(MadaDB_VsearchClusters,
                                         MadaDBdata_filt,
                                         by = "Isolate.ID") 


write_tsv(MadaDB_VsearchClusters_wdata,"./Computational_Resources/Datasets/CultureData/MadaDB_VsearchClusters_wMadaDBdata.txt")
