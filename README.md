# analyzing_peptides

This repository contains scripts that should be used in conjunction with the MultiPep stand-alone program [GitHub link](https://github.com/scheelelab/MultiPep).
We here show a pipeline for using MultiPep predictions conjoined with protein information from UniProt to identify interesting peptides with certain bioactivities from peptidomics data. As an example, we here compare a brain peptidomics dataset with a plasma peptidomics dataset. 


The repository contains two executable scripts:
1.	`getting_peptides.py`: retrieves and cleans the peptides from the mzID files, readies them for MultiPep predicting and print them to TXT files (see Note 1).
2.	`process_predictions.py`: handles MultiPep predicted peptides, adds additional filtering, downloads UniProt protein information and writes output Excel tables and plots distributions of bioactive peptides found in the applied datasets.

The `getting_peptides.py` script utilize a list of known contaminants (“contaminants.fasta”) from [MaxQuant](http://www.coxdocs.org/doku.php?id=maxquant:start_downloads.htm) downloaded 05/04/2022 to remove contaminants. The `process_peptides.py` script uses the MultiPep training data available in a pickle file (`total_classes_seq32.pkl`), to clean the data.

### Retrieving peptides from mzID files
To retrieve peptides from mzID files, one can use the following command:
- `python getting_peptides.py -input_folders brain_PXD008795 plasma_PXD003533 -cont_path contaminants.fasta`

`brain_PXD008795` and `plasma_PXD003533` are the folders containing the mzID files.
In the peptide-retrieval process, peptides will be filtered away based on the following properties:
1.	did not pass the threshold defined when the search algorithm (Mascot or Peaks) was run
2.	peptide was a decoy
3.	peptides contained other symbols than natural amino acids (A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V) 
4.	peptides had a length < 2 or > 200 amino acids.
5.	peptides matched or reversely matched common contamination proteins.
Point 3 will ensure exclusion of peptides with post-translational modifications. Point 3. and 4. were implemented to find peptides that could readily be predicted by MultiPep. 

The script will output `brain_PXD008795.txt`, `plasma_PXD003533.txt` and `brain_PXD008795__plasma_PXD003533.txt` with the detected peptides.

### Predicting with MultiPep
The three TXT files with peptides from the brain and plasma peptidomics data can be predicted with the MultiPep stand-alone program using the following commands:

For the brain peptide list:
-	`python MultiPep_predict_xtra_output.py -input_file brain_PXD008795.txt -output_file brain_PXD008795.csv`

For the plasma peptide list:
-	`python MultiPep_predict_xtra_output.py -input_file plasma_PXD003533.txt -output_file plasma_PXD003533.csv`
For the intersection list:
-	`python MultiPep_predict_xtra_output.py -input_file brain_PXD008795__plasma_PXD003533.txt -output_file brain_PXD008795__plasma_PXD003533.csv`

### Finding high-scoring neuropeptides and peptide hormones
Images of the datasets' predicted bioactivity distributions and tables with the top five highest ranking neurpeptides and peptide hormones can be found using the following command:

-	`python process_predictions.py -input_unique brain_PXD008795.csv plasma_PXD003533.csv -input_inter brain_PXD008795__plasma_PXD003533.csv -multipep_seqs total_classes_seq32.pkl`

The script works via a few steps:
1.	The CSV output files from MultiPep are read and peptides that exist in the MultiPep training data are removed.
2.	The number of peptides with at least one prediction score above 0.5, 0.7 and 0.9, respectively, are plotted for the different datasets.
3.	UniProt protein information is downloaded, and peptides are mapped to protein sequences to find the protein(s) from which they originate.
4.	Tables with top ranking neuropeptides and peptide hormones conjoined with UniProt information are written to Excel files. These tables are additionally cleaned for predicted toxic and hemolytic peptides. 

Peptides that could be found in the MultiPep training data are removed as an attempt to only find peptides not previously described elsewhere. 
