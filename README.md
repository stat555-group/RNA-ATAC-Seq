# Final Project for STAT555

cell lines (1=HSC, 2=CMP, 3=CFUE, 4=Erythroblast)

Assigned Sample pairs:
24; 12; 23; 14; 34; 13;

# link to report
https://docs.google.com/document/d/1E3Oe7eAx9Ex2GmVUrkkG-beus-a870KMqWOrhu4g0-A/edit?usp=sharing

## Downloading the Data 
First, create a ```data``` directory that contains the folders ```rna_seq``` and ``atac_seq``.

```
mkdir data && mkdir data/rna_seq && mkdir data/atac_seq
```
The RNA-Seq data were downloaded in tsv format using the following code:

```
pip install wget
python3 download_data.py -i encode_ids.txt -o data/rna_seq/ -t rna-seq
```

The ATAC-Seq data were downloaded in bigBed format using the following code:

```
python3 download_data.py -i encode_ids.txt -o data/atac_seq/ -t atac-seq
```

For conversion of the ATAC-Seq to tab files, the following scrip was run
while inside the ```atac-seq``` directory:
```
bash convert_atac_seq.sh ../atac_seq/ ../data/atac_seq/
```
