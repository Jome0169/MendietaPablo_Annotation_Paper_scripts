#Wed Jan 22 10:39:57 EST 2020
#Generate to replce pull_regions_assemble_regions.py. Was having issues with
#this earlier where reads from different chromosomes were somehow creaping into
#the output BAM files. This was causing a mayor headache, and the solution below
#feels like a more fool proof soluttion than what I'm doing

### EXAMPLE COMMAND: bash scripts/pull_regions_assemble_regions.sh 00.data/01.bed_files/novel_genes.bed bam_file_list.txt 00.data/03.Isoalted_regions_BAM/novel_regions novel_region_list.txt


bed_file=${1}
bam_file_list=${2}
output_dir=${3}
output_config_file=${4}

echo "region_name file" > ${4}

while read chrom start stop name; do 

    counter=0
    #1:1042000-1042010
    generate_coordinate_query=${chrom}":"${start}"-"${stop}
    generate_coordinate_file_name=${chrom}"_"${start}"-"${stop}".sam"
    generate_coordinate_file_name_bam=${chrom}"_"${start}"-"${stop}".bam"

    echo "Working on Coordinates" ${generate_coordinate_query}

    for i in $(cat ${bam_file_list}) ; do

        if [ ${counter} -eq 0 ]
        then
            #Do the first intersection and append the header
            samtools view -@ 5 -h ${i} ${generate_coordinate_query} > ${3}/${generate_coordinate_file_name}
            #Increase counter by one so no additioanl header
            counter+=1 
        else

            #Don't add the header, but append the reads to the sam file
            samtools view -@ 5 ${i} ${generate_coordinate_query} >> ${3}/${generate_coordinate_file_name}
        fi

    done 

    samtools view -S -b ${3}/${generate_coordinate_file_name} > ${3}/${generate_coordinate_file_name_bam}
    
    #Clean up and remove the SAM file

    rm ${3}/${generate_coordinate_file_name}

    #Append region and BAM file to config file
    echo ${chrom}"_"${start}"-"${stop} ${4}/${generate_coordinate_file_name_bam} >> ${4}

done < ${bed_file} 



