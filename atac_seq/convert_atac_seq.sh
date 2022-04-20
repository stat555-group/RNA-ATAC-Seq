###############################################################################
# This script will process the bigBed files to get the tab files needed
# for downstream analysis.
################################################################################

# Authors: Kobie Kirven, Avantika Diwadkar 


# Download the approprite scripts for transforming the data
wget http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/bigBedToBed
wget http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/bigWigAverageOverBed
wget http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/bedGraphToBigWig

# Download the chromosome sizes file
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes

# Make the downloaded scripts executable
chmod +x bigBedToBed
chmod +x bigWigAverageOverBed
chmod +x bedGraphToBigWig

wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes

# Get the path to the bigBed files
bigBed_path=$2

files=$(ls ${bigBed_path}*.bigBed)

for file in ${files}; do
    # Get the name of the file
    name=$(basename ${file})

    # Get the name of the file without the extension
    name_no_ext=$(echo ${name} | cut -d. -f1)

    # Get the name of the file without the extension and the path
    name_no_ext_no_path=$(echo ${name_no_ext} | cut -d/ -f2)

    # Get the path to the bed file
    bed_path=$(echo ${name_no_ext_no_path} | cut -d. -f1)

    # Convert the bigBed file to a bed file
    ./bigBedToBed ${file} ${bed_path}.bed

    # Convert bed file to bedgraph file.
    cut -f1,2,3,4 ${bed_path}.bed | sort -k1,1 -k2,2n > ${bed_path}.sorted.bedGraph

done


for file in ${files}; do
    # Get the name of the file
    name=$(basename ${file})

    # Get the name of the file without the extension
    name_no_ext=$(echo ${name} | cut -d. -f1)

    # Get the name of the file without the extension and the path
    name_no_ext_no_path=$(echo ${name_no_ext} | cut -d/ -f2)

    # Get the path to the bed file
    bed_path=$(echo ${name_no_ext_no_path} | cut -d. -f1)
    # Solve the problem 
    cat ${bed_path}.sorted.bedGraph | awk -F '\t' -v OFS='\t' '{print $1,$2,$3-1,$4}' > ${bed_path}.sorted.bedGraph.tmp1
    bedtools merge -i ${bed_path}.sorted.bedGraph.tmp1 -c 4 -o mean > ${bed_path}.sorted.bedGraph.tmp1.M.tmp1
    cat ${bed_path}.sorted.bedGraph.tmp1.M.tmp1 | awk -F '\t' -v OFS='\t' '{print $1,$2,$3+1,$4}' > ${bed_path}.sorted.bedGraph

    # Convert bedGraph to bigWig.
    ./bedGraphToBigWig ${bed_path}.sorted.bedGraph mm10.chrom.sizes ${bed_path}.sorted.bw

done


bigBed_path=$1

### use all of the interested peak files to get the master_peak_list. Replace the pk1.bed, pk2.bed ... by the bed files generated from bedgraph files
cat *.bed | sort -k1,1 -k2,2n > ${bigBed_path}all_the_bed_files.pooled.bed
bedtools merge -i ${bigBed_path}all_the_bed_files.pooled.bed > ${bigBed_path}all_the_bed_files.merged.bed
cat ${bigBed_path}all_the_bed_files.merged.bed | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$1"_"$2"_"$3}' > ${bigBed_path}master_peak_list.withpkname.bed


# Get the path to folder with sorted bigwig files and master peak list bed file
bw_path=$1 

files=$(ls ${bw_path}*.bw)

for file in ${files}; do
    # Get the name of the file
    f=$(basename ${file})
    name=$(echo ${f} | cut -d. -f1)
	
    ./bigWigAverageOverBed ${name}.sorted.bw  ${bw_path}master_peak_list.withpkname.bed ${name}.master_peak_list.bigbedfile1.tab
    
    #(Noted: the row will be shuffled so you need to reorder the row)
    sort -k1,1 ${name}.master_peak_list.bigbedfile1.tab > ${name}.master_peak_list.bigbedfile1.pkidsort.tab

done

