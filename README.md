# Protein_coverage_MS description
Calculation of theoretically possible coverage by mass spectrometry for proteins.

# Background
This script was developped as part of the OpenProt pipeline (www.openprot.org). OpenProt is the first proteogenomic resource supporting a polycistronic annotation model for eukaryotic genomes. It provides a deeper annotation of open reading frames (ORFs) while mining experimental data for supporting evidence using cutting-edge algorithms (Brunet MA et al., Nucleic Acids Research, 2019). 
OpenProt predicts all possible ORFs within the transcriptome retrieved from two annotations (NCBI RefSeq and Ensembl) and then retrieves supporting evidence. The prediction pipeline does not enforce a maximal length threshold, although it does filter for a minimal length of 30 codons and an AUG initiating codon. The predicted encoded proteins are then categorized as follow: RefProt or reference proteins are known proteins annotated in NCBI RefSeq, Ensembl and/or UniProt; Novel Isoforms are unannotated proteins with a significant sequence identity to a RefProt from the same gene; and AltProts are unannotated proteins with no significant identity to a RefProt from the same gene. Finally, to assert confidence of the predicted proteins, OpenProt retrieves evidence for each annotated proteins. These are in silico evidence (conservation and prediction of functional domains) and experimental evidence (translation evidence from ribosome profiling data and expression evidence from mass spectrometry data). OpenProt thus offers a deep annotation of the genome of 10 species, identifying novel isoforms and novel proteins in an unbiased yet data-driven manner.
For more information, see the OpenProt publications:
- Brunet MA et al., Nucleic Acids Research, 2019 : https://academic.oup.com/nar/article/47/D1/D403/5123790
- Brunet MA & Roucou X, Journal of Visualized Experiments, 2019 : https://www.jove.com/video/59589/mass-spectrometry-based-proteomics-analyses-using-openprot-database
- Brunet MA et al., Current Protocols in Bioinformatics, 2020 : https://currentprotocols.onlinelibrary.wiley.com/doi/full/10.1002/cpbi.103

# Details of calculations
To help the OpenProt users with the interpretation of a MS score of 0, we added MS coverage statistics to the MS tab of the details page for each protein. Amongst these metrics is the theoretical coverage for each given protein. The theoretical coverage is calculated from all possible tryptic peptides that would fit the OpenProt criteria to be assigned: a minimal length of 7 amino acids, a maximal mass of 4,600 Da and peptide unicity given the protein type (RefProt, Novel Isoform or AltProt).

