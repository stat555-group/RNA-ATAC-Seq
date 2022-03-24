# Final Project for STAT555

cell lines (1=HSC, 2=CMP, 3=CFUE, 4=Erythroblast)

Assigned Sample pairs:
24; 12; 23; 14; 34; 13; 


## Downloading the Data 
First, create a ```data``` directory that contains the folders ```rna_seq``` and ``atac_seq``.

```
mkdir data && mkdir data/rna_seq && mkdir atac_seq
```
The RNA-Seq data were downloaded in tsv format using the following code:

```
python3 download_data.py -i encode_ids.txt -o data/rna_seq/ -t rna-seq
```

The ATAC-Seq data were downloaded in bigBed format using the following code:

```
python3 download_data.py -i encode_ids.txt -o data/atac_seq/ -t atac-seq
```