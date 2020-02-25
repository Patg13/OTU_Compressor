# OTU_Compressor
Program to compress spurious OTUs using GenBank Accession number 

Besthit Blast Filter & Compressor V1
By: Patrick Gagne

--------------------------------------------

NOTICE:
This program regroup besthits lines using Accession Number ONLY
It doesn't consider % Identity or BitScore

The Alignment length is used for filtration purpose
and HAVE PRIORITY over regrouping

-------------------------------------------

usage: Besthits_CNF.py [-h] -in BESTHIT_FILE [-aln_thr ATHRES_VALUE]
                       [-report REPORT_FILENAME] [-checkall] -out
                       OUTPUT_FILENAME

Besthit Blast Filter & Compressor V1

optional arguments:
  -h, --help            show this help message and exit
  
  -in BESTHIT_FILE      Best Hits Blast file using this format (OTU_id Size
                        Iden Acc Piden Aln_Len Bitscore Frequencies_for_each_samples)
                        [REQUIRED]
                  
  -aln_thr ATHRES_VALUE
                        Minimum Alignment Length for hit [default=100]
                        
  -report REPORT_FILENAME
                        Program Report filename (Contain informations about
                        rejected OTUs) [OPTIONAL]
                        
  -checkall             Check validity for the whole file (This will take more
                        time)
                        
  -out OUTPUT_FILENAME  Output Filename [REQUIRED]
