from __future__ import division
from optparse import OptionParser
from Bio import SeqIO
import sys
import subprocess
import os
import cPickle as Pickle
from collections import OrderedDict
import numpy as np
import scipy.stats as stat

CODON_TYPES = ['UUU', 'UUC', 'UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG', 'AUU', 'AUC', 'AUA', 'AUG', 'GUU', 'GUC', 'GUA',
               'GUG', 'UCU', 'UCC', 'UCA', 'UCG', 'CCU', 'CCC', 'CCA', 'CCG', 'ACU', 'ACC', 'ACA', 'ACG', 'GCU', 'GCC',
               'GCA', 'GCG', 'UAU', 'UAC', 'CAU', 'CAC', 'CAA', 'CAG', 'AAU', 'AAC', 'AAA', 'AAG', 'GAU', 'GAC', 'GAA',
               'GAG', 'UGU', 'UGC', 'UGG', 'CGU', 'CGC', 'CGA', 'CGG', 'AGU', 'AGC', 'AGA', 'AGG', 'GGU', 'GGC', 'GGA',
               'GGG', 'UAA', 'UAG', 'UGA']

genetic_code = {'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C', 'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',
                'UUA': 'L', 'UCA': 'S', 'UAA': '*', 'UGA': '*', 'UUG': 'L', 'UCG': 'S', 'UAG': '*', 'UGG': 'W',
                'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R', 'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
                'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R', 'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',
                'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S', 'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
                'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R', 'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',
                'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G', 'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
                'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G', 'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'}

AMINO_ACIDS = ['A', 'R', 'D', 'N', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '*']

MULTIPLE_MAPPED_GENE_READ_COUNTS_FILE = "Multiple_mapped_gene_read_counts.tab"

def processSamFile(folder,
                   bam_file):
    # create the folder
    sam_file = os.path.join(folder, 'samfile.sam')
    
    try:
        subprocess.call('samtools view ' + bam_file + ' > ' + samfile, shell=True)
        
    except Exception:
        print('samtools failed. Check the input bam file.')
        
    return sam_file

def processAnnotationFile(folder,
                         annotation_file):
    file_name, ext = os.path.splittext(annotation_file)
    if ext == 'gff' or ext == 'gff3':
        if 'sacCer3' in annotation_file or 'saccer3' in annotation_file:
            annotation_file = process_gff_sacCer3(annotation_file)
        elif 'ecoli' in annotation_file:
            annotation_file = process_gff_ecoli(annotation_file)
        else:
            print('[Warning]: GFF file entered as input for gene annotations. Exercise caution. Parsing of GFF file is optimized only for sacCer3 and E.coli. Check the format')
            
        return annotation_file
    

def samparser_genome(sfile, frag_min, frag_max, three_prime):
    # Parse the SAM file and quantify mapped reads to chromosome positions.
    dict_count = {}  # This dictionary will contain the count values and other attributes for every genomic position
    dict_mul_count = {}  # This dictionary will contain the count values of multiple mapped reads
    dict_total = {}
    # Multiple mapped reads are counted separately and if they are greater than certain % of total reads mapped on a gene, that gene will not be considered for analysis
    log_file = open('sam_parser_logfile.log', 'w')
    frag_range = frag_max - frag_min + 1
    read_count = 0
    mul_read_count = 0
    if three_prime:
        print('Mapping reads by 3\' end...')
    else:
        print('Mapping reads by 5\' end...')
    counter = 0
    # Open the SAM file and parse the mapped reads according to chromosome positions, read length and strand
    with open(sfile) as sam_file:
        for line in sam_file:
            # Ignore header files
            if line[0] != '@':
                fields = line.split('\t')
                if len(fields) > 19:
                    multi_flag = fields[19].strip()  # Sometimes there are "NH:i:1\n"
                    # If the read has more than one alignment then report it as multiple mapping
                    if multi_flag != 'NH:i:1':
                        multi_align = True
                    else:
                        multi_align = False
                else:
                    multi_align = False
                sam_flag = int(fields[1])
                chr_no = fields[2]
                read_length = len(fields[9])

                if sam_flag in (0, 256):
                    # Assuming no insertions and deletions
                    # If mapping by 3' end
                    if three_prime:
                        position = int(fields[3]) + read_length - 1
                    # If mapping by 5' end
                    else:
                        position = int(fields[3])

                elif sam_flag in (16, 272):
                    # If mapping by 3' end
                    if three_prime:
                        position = int(fields[3])
                    # If mapping by 5' end
                    else:
                        position = int(fields[3]) + read_length - 1

                # Initialize the dict if not done already
                if chr_no not in dict_count:
                    dict_count[chr_no] = {}
                    dict_total[chr_no] = 0
                # Initialize an array of zeroes corresponding to reads of different lengths. Two values for each length representing the positive and negative strands
                # The last two zeroes will correspond to multiple aligned reads according to strand
                if position not in dict_count[chr_no]:
                    dict_count[chr_no][position] = [0] * (frag_range * 2 + 2)

                # If read is multiple mapped, then initialize the multiple count dict as well
                if multi_align and chr_no not in dict_mul_count:
                    dict_mul_count[chr_no] = {}
                if multi_align and position not in dict_mul_count[chr_no]:
                    dict_mul_count[chr_no][position] = [0] * (frag_range * 2 + 2)

                pos_in_value = (read_length - frag_min) * 2
                # Count the read according to its length and strand
                if frag_min <= read_length <= frag_max:
                    try:
                        if sam_flag == 0:  # Primary alignment on forward strand
                            if not multi_align:
                                dict_count[chr_no][position][pos_in_value] += 1
                                dict_total[chr_no] += 1
                                read_count += 1
                            else:
                                # Multiple mapped reads are counted separately
                                dict_count[chr_no][position][-2] += 1
                                dict_mul_count[chr_no][position][pos_in_value] += 1
                                mul_read_count += 1
                        if sam_flag == 16:  # Primary alignment on reverse strand
                            if not multi_align:
                                dict_count[chr_no][position][pos_in_value + 1] += 1
                                dict_total[chr_no] += 1
                                read_count += 1
                            else:
                                # Multiple mapped reads are counted separately. Last two columns are initialized for mul mapped reads for +ve and -ve strands respectively
                                dict_count[chr_no][position][-1] += 1
                                dict_mul_count[chr_no][position][pos_in_value + 1] += 1
                                mul_read_count += 1
                        # Not primary alignment. It will counted under multiple aligned reads
                        if sam_flag == 256:
                            position = int(fields[3])
                            dict_count[chr_no][position][-2] += 1
                            dict_mul_count[chr_no][position][pos_in_value] += 1

                        if sam_flag == 272:  # Not primary alignment and on reverse strand
                            position = int(fields[3]) + read_length - 1
                            dict_count[chr_no][position][-1] += 1
                            dict_mul_count[chr_no][position][pos_in_value + 1] += 1
                    except KeyError:
                        log_file.write("The KeyError line is \n " + line + "\n")
            counter += 1
            sys.stdout.write("Line of SAM file currently being parsed: {0}.\t\r".format(counter))
            sys.stdout.flush()

    print('\nSAM file parsed for total ' + str(read_count) + ' reads mapping onto ' + str(len(dict_count)) + ' chromosomes.'
    print(str(mul_read_count) + ' reads are multiple aligned mapped to ' + str(len(dict_mul_count)) + ' chromosomes.'))

    return dict_count, dict_mul_count


def samparser_transcriptome(sfile, frag_min, frag_max, three_prime):
    # Parse the SAM file and quantify mapped reads to chromosome positions.
    dict_count = {}  # This dictionary will contain the count values and other attributes for every genomic position
    dict_mul_count = {}  # This dictionary will contain the count values of multiple mapped reads
    total_count = {}
    # Multiple mapped reads are counted separately and if they are greater than 0.1% of total reads mapped on a gene, that gene will not be considered for analysis
    log_file = open('sam_parser_logfile.log', 'w')
    frag_range = frag_max - frag_min + 1
    read_count = 0
    mul_read_count = 0
    discarded_reads = 0
    if three_prime:
        print('Mapping reads by 3\' end...')
    else:
        print('Mapping reads by 5\' end...')
    counter = 0
    with open(sfile) as sam_file:
        for line in sam_file:
            # Ignore header files
            if line[0] != '@':
                fields = line.split('\t')
                sam_flag = int(fields[1])
                gene = fields[2]
                # if gene not in filtered_genes:
                #     continue
                if len(fields) > 19:
                    multi_flag = fields[19].strip()  # Sometimes there are "NH:i:1\n"
                    # If the read has more than one alignment then report it as multiple mapping
                    if multi_flag != 'NH:i:1':
                        multi_align = True
                    else:
                        multi_align = False
                else:
                    multi_align = False

                read_length = len(fields[9])
                pos_in_value = (read_length - frag_min)
                if gene not in dict_count:
                    dict_count[gene] = {}
                    total_count[gene] = 0
                # If mapping by 3' end
                if three_prime:
                    position = int(fields[3]) + read_length - 1
                # If mapping by 5' end
                else:
                    position = int(fields[3])
                if position not in dict_count[gene]:
                    dict_count[gene][position] = [0] * (frag_range + 1)
                if multi_align:
                    if gene not in dict_mul_count:
                        dict_mul_count[gene] = {}
                    if position not in dict_mul_count[gene]:
                        dict_mul_count[gene][position] = [0] * (frag_range + 1)

                if frag_min <= read_length <= frag_max:
                    try:
                        if sam_flag == 0:  # Primary alignment on forward strand
                            if not multi_align:
                                dict_count[gene][position][pos_in_value] += 1
                                read_count += 1
                                total_count[gene] += 1
                            else:
                                # Multiple mapped reads are counted separately
                                dict_count[gene][position][-1] += 1
                                dict_mul_count[gene][position][pos_in_value] += 1
                                mul_read_count += 1
                        # Not primary alignment. It will counted under multiple aligned reads
                        elif sam_flag == 256:
                            dict_count[gene][position][-1] += 1
                            dict_mul_count[gene][position][pos_in_value] += 1
                        elif sam_flag == 16:
                            discarded_reads += 1

                    except KeyError:
                        log_file.write("The KeyError line is \n " + line + "\n")
            counter += 1
            sys.stdout.write("Line of SAM file currently being parsed: {0}.\t\r".format(counter))
            sys.stdout.flush()

    print('\nSAM file parsed for total ' + str(read_count) + ' reads mapping onto ' + str(len(dict_count)) + ' genes.')
    print(str(mul_read_count) + ' reads are multiple aligned mapped to ' + str(len(dict_mul_count)) + ' genes.')
    print(str(discarded_reads) + ' reads were mapped spuriously.')
    return dict_count, dict_mul_count, total_count


def create_cds_counts_transcriptome(idx_file, seq_file, output, sam_count_dict, dict_mul_count, total_count, frag_min, frag_max, three_prime, fast_mode=True):

    logfile = open(output+'makecdstable.log', 'w')
    frag_range = frag_max - frag_min + 1
    mul_gene_list = []
    dict_mul_genes = {}
    dict_len = {}
    idx_dict = {}
    dict_count_len = {}
    dict_mul_count_len = {}
    # These dicts will contain dicts of genes and their read counts (as lists) according to fragment length
    for fsize in range(frag_min, frag_max + 1):
        dict_count_len[fsize] = {}

    with open(idx_file) as f:
        for lines in f:
            fields = lines.strip().split('\t')
            # fields = [gene name, CDS start position, CDS length]
            idx_dict[fields[0]] = [int(fields[1]), int(fields[2])]

    seq_dict = {}

    for seq_record in SeqIO.parse(seq_file, "fasta"):
        gene_id = str(seq_record.id).split(' ')[0]
        seq_dict[gene_id] = seq_record.seq

    dict_start = {}
    print('Starting to make the CDS table')
    for gene in sam_count_dict:
        try:
            start_pos, cds_len = idx_dict[gene]
        except KeyError:
            logfile.write('No index available for gene ' + gene + '\n')
            continue
        # gene_count_file.write(gene + '\t' + str(total_count[gene]) + '\t' + str(cds_len) + '\n')

        if start_pos in [0, -1]:
            logfile.write('No UTR region for gene ' + gene + '\n')
            continue
        # If a gene has very sparse reads, it is better to leave it out as it will not meet the filtering criteria.
        # below criteria states that avg is less than 1. This is done for faster processing.
        if fast_mode and total_count[gene] < cds_len:
            logfile.write('SPARSELY POPULATED - Being filtered:\t' + gene + '\t' + str(total_count[gene]) + '\t' + str(cds_len) + '\n')
            continue
        dict_len[gene] = cds_len
        multi_genes = 'N'
        start_utr = -start_pos
        cds_pos = start_utr
        for pos in range(1, len(seq_dict[gene]) + 1):
            # for pos in sorted(sam_count_dict[gene]):
            nuc = seq_dict[gene][pos - 1]
            # We are not interested in positions beyond 50 nt around the CDS region. Hence we skip these positions
            if cds_pos < -50 or cds_pos > cds_len + 50:
                cds_pos += 1
                continue
            if cds_pos == 0:
                cds_pos += 1
            try:
                count_list = sam_count_dict[gene][pos][0:frag_range]
            except KeyError:
                count_list = [0] * frag_range

            try:
                multi_aligned = sam_count_dict[gene][pos][-1]
            except KeyError:
                multi_aligned = 0

            try:
                if multi_aligned > 0:
                    mul_count_list = dict_mul_count[gene][pos][0:frag_range]
            except KeyError:
                mul_count_list = [0] * frag_range

            # We will create dictionaries with fragment length as keys and dict of gene with list of read counts at every position as our values
            for fsize in range(frag_min, frag_max + 1):
                if gene not in dict_count_len[fsize]:
                    dict_count_len[fsize][gene] = []
                    dict_start[gene] = cds_pos
                dict_count_len[fsize][gene].append(count_list[fsize - frag_min])

            if multi_aligned > 0 and gene not in mul_gene_list:
                mul_gene_list.append(gene)
                if gene not in dict_mul_genes:
                    dict_mul_genes[gene] = {}
            if multi_aligned > 0:
                dict_mul_genes[gene][cds_pos] = ['X', cds_len, nuc, 'transcriptome', pos, multi_aligned, multi_genes] + mul_count_list
                # As done for actual counts, multiple mapped reads are counted according to fragment length
                for fsize in range(frag_min, frag_max + 1):
                    if gene not in dict_mul_count_len:
                        dict_mul_count_len[gene] = {}
                    if fsize not in dict_mul_count_len[gene]:
                        dict_mul_count_len[gene][fsize] = {0: 0, 1: 0, 2: 0}
                    if cds_pos < 0:
                        frame = cds_pos % 3
                    else:
                        frame = cds_pos % 3 - 1
                        if frame == -1:
                            frame = 2
                    if three_prime:
                        pos_range = range(1, cds_len + frag_max)
                    else:
                        pos_range = range(-frag_max, cds_len)
                    if cds_pos in pos_range:
                        dict_mul_count_len[gene][fsize][frame] += mul_count_list[fsize - frag_min]

            cds_pos += 1

    # Write out the number of mutliple mapped reads for each gene according to fragment size and frame. For each fragment size, write the reads for frame 0, 1 and 2.
    mul_out_file = open(output+MULTIPLE_MAPPED_GENE_READ_COUNTS_FILE, 'w')

    for gene in dict_mul_count_len:
        mul_out_file.write(str(gene))
        for fsize in range(frag_min, frag_max + 1):
            for frame in xrange(3):
                mul_out_file.write('\t' + str(dict_mul_count_len[gene][fsize][frame]))
        mul_out_file.write('\n')
    mul_out_file.close()

    for fsize in range(frag_min, frag_max + 1):
        count_file = open(output+'Read_counts_' + str(fsize) + '.tab', 'w')
        for gene, reads_list in dict_count_len[fsize].iteritems():
            start_pos, cds_len = idx_dict[gene]
            start_idx = dict_start[gene]
            # Writing the read count for transcripts
            count_file.write(gene + '\t' + str(start_idx) + '\t' + str(cds_len) + '\t' + ','.join(map(str, reads_list)) + '\n')
        count_file.close()


def create_cds_counts_genome(annotation_file, genome, output, sam_parsed_count_dict, dict_mul_count, frag_min, frag_max, three_prime, extra_overlap, remove_overlapped_genes=True):
    # This module will take the count dictionary parsed from the sam file containing the read counts at each genomic position and map it to gene positions

    complimentarydict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}

    dict_gene, dict_cds_count, dict_cds_info, genome_dict, overlap_genes, dict_len = cdsparser(annotation_file, genome, extra_overlap=extra_overlap)

    dict_count_len = {}
    dict_mul_count_len = {}
    # These dicts will contain dicts of genes and their read counts (as lists) according to fragment length
    for fsize in range(frag_min, frag_max + 1):
        dict_count_len[fsize] = {}

    counter = 0

    # For every gene in the dataset
    for gene_name in dict_cds_count:
        chr_num, strand = dict_gene[gene_name]
        frag_range = frag_max - frag_min + 1
        try:
            gene_length = dict_len[gene_name]
        except KeyError:
            print(dict_cds_info[gene_name])
            print(dict_len[gene_name])
            print('KeyError in length calculation for dict_cds_info in gene ' + str(gene_name))

        # dict_cds_info contains lists of CDS start and end positions as a list. For example if there are two exons for a gene X
        # dict_cds_info['X'] = [[111234, 111345],[112122,112543]]
        # We get reads at 50 nt around the first cds position and the last cds position and update these lists
        dict_cds_info[gene_name][0][0] -= 50
        dict_cds_info[gene_name][-1][1] += 50
        # Gene position index
        gene_position_start = -50
        # For -ve strand, this will be start of position as they go in opposite direction
        gene_position_end = gene_length + 50
        for k in dict_cds_info[gene_name]:
            for pos in range(k[0], k[1] + 1):
                if gene_position_start == 0:
                    gene_position_start += 1
                if gene_position_end == 0:
                    gene_position_end += -1
                try:
                    nuc = genome_dict[chr_num][pos - 1]
                except KeyError:
                    print('KEYERROR in genome dict for chromosome ' + str(chr_num))
                    print('Genome dict is ' + str(genome_dict.keys()))
                except IndexError:
                    # A gene maybe present at the extreme end of chromosome (in E.coli) and hence 50 positions downstream will not be present
                    nuc = 'N'

                try:
                    counts_list = sam_parsed_count_dict[chr_num][pos][0:frag_range * 2]
                except KeyError:
                    # sam_parsed_count_dict is initialized only for genes with read counts.
                    # Positions with zero reads will lead to KeyErrors and we assign zero value to this position
                    counts_list = [0] * (frag_range * 2)
                # check if this position has multiple aligned reads
                if chr_num not in dict_mul_count:
                    pass
                elif pos in dict_mul_count[chr_num]:
                    mul_count_list = dict_mul_count[chr_num][pos][0:frag_range * 2]

                if strand == '+':
                    try:
                        multi_mapped = sam_parsed_count_dict[chr_num][pos][-2]
                    except KeyError:
                        multi_mapped = 0
                    # Get the read counts on the positive strand
                    strand_counts = countstrand(counts_list, '+')
                    if multi_mapped > 0:
                        # If there are multiple mapped reads at this position, get the read counts for those multiple mapped reads
                        mul_strand_counts = countstrand(mul_count_list, '+')
                    # For positive strand, the current position in "for" loop is from gene start side
                    current_pos = gene_position_start

                if strand == '-':
                    try:
                        multi_mapped = sam_parsed_count_dict[chr_num][pos][-1]
                    except KeyError:
                        multi_mapped = 0

                    # Get the read counts on the negative strand
                    strand_counts = countstrand(counts_list, '-')
                    if multi_mapped > 0:
                        mul_strand_counts = countstrand(mul_count_list, '-')
                    # Since the strand is negative, we take the complement of the current nucleotide
                    nuc = complimentarydict[nuc]
                    current_pos = gene_position_end

                # We will create dictionaries with fragment length as keys and dict of gene with list of read counts at every position as our values
                for fsize in range(frag_min, frag_max + 1):
                    if gene_name not in dict_count_len[fsize]:
                        dict_count_len[fsize][gene_name] = []
                    dict_count_len[fsize][gene_name].append(strand_counts[fsize - frag_min])

                if multi_mapped > 0:
                    mul_mapped_reads = 0
                    # As done for actual counts, multiple mapped reads are counted according to fragment length
                    for fsize in range(frag_min, frag_max + 1):
                        if gene_name not in dict_mul_count_len:
                            dict_mul_count_len[gene_name] = {}
                        if fsize not in dict_mul_count_len[gene_name]:
                            dict_mul_count_len[gene_name][fsize] = {0: 0, 1: 0, 2: 0}
                        if current_pos < 0:
                            frame = current_pos % 3
                        else:
                            frame = current_pos % 3 - 1
                            if frame == -1:
                                frame = 2
                        if three_prime:
                            pos_range = range(1, gene_length + frag_max)
                        else:
                            pos_range = range(-frag_max, gene_length)
                        if current_pos in pos_range:
                            dict_mul_count_len[gene_name][fsize][frame] += mul_strand_counts[fsize - frag_min]
                            mul_mapped_reads += mul_strand_counts[fsize - frag_min]
                    if mul_mapped_reads == 0:
                        del dict_mul_count_len[gene_name]

                if strand == '-':
                    gene_position_end += -1
                if strand == '+':
                    gene_position_start += 1
        counter += 1
        sys.stdout.write("Out of " + str(len(dict_cds_count)) + " transcripts, currently processing transcript {0}.\t\r".format(counter))
        sys.stdout.flush()

    # Determine the percentage of reads which are multiple mapped to a gene and discard it if it is greater than the threshold set for multiple map filter.
    # This is done specific to each fragment size
    mul_out_file = open(output+MULTIPLE_MAPPED_GENE_READ_COUNTS_FILE, 'w')

    for gene in dict_mul_count_len:
        mul_out_file.write(gene)
        for fsize in range(frag_min, frag_max + 1):
            for frame in xrange(3):
                mul_out_file.write('\t' + str(dict_mul_count_len[gene][fsize][frame]))
        mul_out_file.write('\n')
    mul_out_file.close()

    for fsize in range(frag_min, frag_max + 1):
        count_file = open(output+'Read_counts_' + str(fsize) + '.tab', 'w')
        for gene, reads_list in dict_count_len[fsize].iteritems():
            chr_num, strand = dict_gene[gene]
            length = dict_len[gene]
            if remove_overlapped_genes and gene in overlap_genes:
                continue
            if strand == '-':
                # For negative strand, read count were appended from the opposite end
                reads_list = reversed(reads_list)
            count_file.write(gene + '\t-50\t' + str(length) + '\t' + ','.join(map(str, reads_list)) + '\n')
        count_file.close()


# Calculate total length of multiCDS region
def multicdslen(cds_positions):
    length_total = 0
    for each_pos in cds_positions:
        each_len = each_pos[1] - each_pos[0] + 1
        length_total += each_len
    return length_total


# Count the reads on 5' end with same strand
def countstrand(counts_list, strand):
    frag_range = int(len(counts_list) / 2)
    strand_counts = []
    # initialize the strand count list
    for q in range(0, frag_range):
        strand_counts.append(0)

    # Start filling the strand count list with actual values from the count_list
    for r in range(0, frag_range):
        if strand == '+':
            strand_counts[r] = counts_list[r * 2]
        if strand == '-':
            strand_counts[r] = counts_list[r * 2 + 1]
    return strand_counts


# Parse the genome reference file
def read_fasta(genome_file):
    chrom, sequence = None, []
    for line in genome_file:
        line = line.rstrip()  # Removing any trailing whitespace characters
        if line.startswith(">"):
            if chrom:
                yield (chrom, ''.join(sequence))
            chrom, sequence = line, []
        else:
            sequence.append(line)  # Appends each sequence line to list. Finally joins in at the end of chromosome
    if chrom:
        yield (chrom, ''.join(sequence))


# Parse the GFF file to get the positions and no. of CDS regions in a gene.
# Important data structures: dictcdscount and dictcdsinfo
def process_gff_sacCer3(gff):
    dictgene = {}  # Contains gene names as keys and their start and end positions as a list of values
    dictcdscount = {}  # Contains gene names as keys and no. of cds regions(exons) for that particular gene
    dictcdsinfo = {}  # Contains the start and end positions of the cds regions as list of lists
    dict_len = {}

    with open(gff) as gff_file:
        for line in gff_file:
            line_list = line.strip().split('\t')
            if line.startswith('chr'):
                left_pos = int(line_list[3])
                right_pos = int(line_list[4])
                chrnum = line_list[0]
                strand = line_list[6]
                if line_list[2] == 'CDS':
                    cds_annotation = line_list[8]
                    # This is for sacCer3 genome annotation file
                    cds_name = cds_annotation.split(';')[0].split('=')[1].split('_')[0]
                    # cds_name = cds_annotation.split(';')[0].split('=')[1]
                    if cds_name not in dictcdscount:
                        dictcdscount[cds_name] = 1
                        dictcdsinfo[cds_name] = []
                    else:
                        dictcdscount[cds_name] += 1
                    dictgene[cds_name] = [chrnum, strand]
                    dictcdsinfo[cds_name].append([left_pos, right_pos])
    # intron_genes_file = open("list_genes_introns.tab", 'w')
    cds_file = open("cds_info.tab", "w")

    for gene in dictcdsinfo:
        chrnum, strand = dictgene[gene]
        dict_len[gene] = multicdslen(dictcdsinfo[gene])
        #  Gene name  Chr number Strand  Length of gene (CDS regions only) No of CDS regions
        cds_file.write(gene + '\t' + chrnum + '\t' + strand + '\t' + str(dict_len[gene]) + '\t' + str(dictcdscount[gene]))
        # Start and end of each CDS region
        for cds in dictcdsinfo[gene]:
            cds_file.write('\t' + str(cds[0]) + '\t' + str(cds[1]))
        cds_file.write('\n')

    # for gene in dictcdscount:
    #     if dictcdscount[gene] > 1:
    #         intron_genes_file.write(gene + '\n')
    gff_file.close()
    cds_file.close()
    return cds_file.name


def process_gff_ecoli(gff):
    dictgene = {}  # Contains gene names as keys and their start and end positions as a list of values
    dictcdsinfo = {}  # Contains the start and end positions of the cds regions as list of lists
    dict_psuedo = {}
    dict_len = {}
    problem_genes = []
    chrom = ''
    with open(gff) as gff_file:
        for line in gff_file:
            line_list = line.strip().split('\t')
            if not line.startswith('##'):
                try:
                    left_pos = int(line_list[3])
                    right_pos = int(line_list[4])
                    strand = line_list[6]
                    if left_pos == 3698003:
                        print(left_pos, right_pos, strand)
                    chrom = line_list[0]
                    if line_list[2] == 'gene':
                        gene_anno = line_list[8]
                        fields = gene_anno.split(';')
                        gene_info = {}
                        for field in fields:
                            fd = field.split('=')
                            gene_info[fd[0]] = fd[1]
                        gene_name = gene_info['ID']
                        if 'pseudo' in gene_info:
                            if gene_info['pseudo'] == 'true':
                                dict_psuedo[gene_name] = [left_pos, right_pos, strand]
                                continue
                            else:
                                dictgene[gene_name] = [left_pos, right_pos, gene_info]
                        else:
                            dictgene[gene_name] = [left_pos, right_pos, gene_info]
                    if line_list[2] == 'CDS':
                        cds_annotation = line_list[8]
                        fields = cds_annotation.split(';')
                        cds_info = {}
                        for field in fields:
                            fd = field.split('=')
                            cds_info[fd[0]] = fd[1]
                        gene_alias = cds_info['Parent']
                        if gene_alias in dict_psuedo or 'pseudo' in cds_info:
                            continue
                        if gene_alias in dictgene:
                            cds_name = dictgene[gene_alias][2]['locus_tag']

                        if gene_alias in dictgene:
                            if left_pos != dictgene[gene_alias][0]:
                                print('CDS for gene ' + cds_name + ' does not start at gene start. The cds start is at ' + str(left_pos) + ' and gene start is at ' + str(
                                    dictgene[gene_alias][0]))
                                print('CDS for gene ' + cds_name + ' does not stop at gene stop. The cds stop is at ' + str(right_pos) + ' and gene stop is at ' + str(
                                    dictgene[gene_alias][1]))
                                problem_genes.append(cds_name)
                        else:
                            print('Gene with only CDS annotation is ' + cds_name)
                        if cds_name not in dictcdsinfo:
                            dictcdsinfo[cds_name] = [left_pos, right_pos, strand]
                        else:
                            print('Second CDS present for gene ' + cds_name)
                            print('First CDS')
                            print(dictcdsinfo[cds_name])
                            print('Second CDS')
                            print(left_pos, right_pos, strand)
                            dictcdsinfo[cds_name] = [left_pos, right_pos, strand]
                        gene_length = right_pos - left_pos + 1
                        dict_len[cds_name] = gene_length
                except IndexError:
                    print('Weird line:' + line)
    outfile = open("../data_files/ecoli/CDS_info.tab", 'w')
    for gene in problem_genes:
        del dictcdsinfo[gene]
    for gene in dictcdsinfo:
        left_pos, right_pos, strand = dictcdsinfo[gene]
        outfile.write(gene + '\t' + str(chrom) + '\t' + str(strand) + '\t' + str(dict_len[gene]) + '\t1\t' + str(left_pos) + '\t' + str(right_pos) + '\n')

    return outfile


def cdsparser(annot_file, genome, extra_overlap=0):
    dict_cds_count = {}
    dict_cds_info = {}
    dict_len = {}
    dict_gene_info = {}

    with open(annot_file) as f:
        for lines in f:
            fields = lines.strip().split('\t')
            gene, chrnum, strand, length, count = fields[:5]
            dict_cds_count[gene] = int(count)
            dict_len[gene] = int(length)
            dict_gene_info[gene] = [chrnum, strand]
            dict_cds_info[gene] = []
            for i in range(5, len(fields), 2):
                dict_cds_info[gene].append([int(fields[i]), int(fields[i + 1])])

    genomedict = {}  # This dictionary contains chromosomes (chrI, chrII, etc) as keys and the sequence as its value
    with open(genome) as fp:
        for name, seq in read_fasta(fp):
            genomedict[name[1::]] = seq
    # Get the list of genes which have overlapping cds regions
    overlap_genes = find_overlapping_genes(dict_cds_info, dict_gene_info, extra_len=extra_overlap)

    return dict_gene_info, dict_cds_count, dict_cds_info, genomedict, overlap_genes, dict_len


# Find overlapping genes from a CDS dictionary of genes
def find_overlapping_genes(dict_cds_info, dict_gene, extra_len=0):
    dict_search_overlap = {}
    dict_overlap = {}
    overlap_list = []

    for gene_name in dict_cds_info:
        chrnum, strand = dict_gene[gene_name]
        # Left most position
        leftpos = int(dict_cds_info[gene_name][0][0])
        # Right most position
        rightpos = int(dict_cds_info[gene_name][-1][1])
        # Comparing the boundaries of individual transcripts
        try:
            regions = dict_search_overlap[chrnum]
            for gene in regions:

                if strand == gene[3]:
                    if leftpos - extra_len < gene[2] and rightpos + extra_len > gene[1]:
                        dict_overlap[gene[0]] = [chrnum, gene[1], gene[2], strand, gene_name, leftpos, rightpos, strand]
                        dict_overlap[gene_name] = [chrnum, leftpos, rightpos, strand, gene[0], gene[1], gene[2], strand]
            dict_search_overlap[chrnum].append([gene_name, leftpos, rightpos, strand])
        except KeyError:
            dict_search_overlap[chrnum] = [[gene_name, leftpos, rightpos, strand]]

    for gene, info in dict_overlap.iteritems():
        overlap_list.append(gene)
    overlap_list.sort()

    return overlap_list



def get_transcript_sequences(annotation_file, genome, output, extra_overlap=0):
    # This module will take the count dictionary parsed from the sam file containing the read counts at each genomic position and map it to gene positions

    complimentarydict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}

    dict_gene, dict_cds_count, dict_cds_info, genome_dict, overlap_genes, dict_len = cdsparser(annotation_file, genome, extra_overlap=extra_overlap)
    print('Parsed the annotation and genome files')
    print('Size of dict_Cds_Count dict is '+str(len(dict_cds_count)))
    print('Size of genome dict'+str(len(genome_dict)))
    counter = 0
    nuc_dict = {}
    # For every gene in the dataset
    for gene_name in dict_cds_count:
        chr_num, strand = dict_gene[gene_name]
        try:
            gene_length = dict_len[gene_name]
        except KeyError:
            print(dict_cds_info[gene_name])
            print(dict_len[gene_name])

        nuc_dict[gene_name] = []
        # dict_cds_info contains lists of CDS start and end positions as a list. For example if there are two exons for a gene X
        # dict_cds_info['X'] = [[111234, 111345],[112122,112543]]
        # We get reads at 50 nt around the first cds position and the last cds position and update these lists
        dict_cds_info[gene_name][0][0] -= 50
        dict_cds_info[gene_name][-1][1] += 50
        # Gene position index
        gene_position_start = -50
        # For -ve strand, this will be start of position as they go in opposite direction
        gene_position_end = gene_length + 50
        for k in dict_cds_info[gene_name]:
            for pos in range(k[0], k[1] + 1):
                if gene_position_start == 0:
                    gene_position_start += 1
                if gene_position_end == 0:
                    gene_position_end += -1
                try:
                    nuc = genome_dict[chr_num][pos - 1]
                except KeyError:
                    print('KEYERROR in genome dict for chromosome ' + str(chr_num))
                    print('Genome dict is ' + str(genome_dict.keys()))
                except IndexError:
                    # A gene maybe present at the extreme end of chromosome (in E.coli) and hence 50 positions downstream will not be present
                    nuc = 'N'

                if strand == '-':
                    # Since the strand is negative, we take the complement of the current nucleotide
                    nuc = complimentarydict[nuc]
                nuc_dict[gene_name].append(nuc)
        print('Got gene sequence for '+gene_name)
        counter += 1
        sys.stdout.write("Out of " + str(len(dict_cds_count)) + " transcripts, currently processing transcript {0}.\t\r".format(counter))
        sys.stdout.flush()

        count_file = open(output+'Transcript_sequence.tab', 'w')
        for gene in nuc_dict:
            chr_num, strand = dict_gene[gene]
            length = dict_len[gene]
            if strand == '+':
                sequence = nuc_dict[gene]
            if strand == '-':
                # For negative strand, read count were appended from the opposite end
                sequence = reversed(nuc_dict[gene])
            count_file.write(gene + '\t-50\t' + str(length) + '\t' + ''.join(map(str, sequence)) + '\n')
        count_file.close()


def generate_asite_profiles(frag_min, frag_max, offfile, infolder):
    offsets = {}
    with open(offfile) as f:
        for lines in f:
            fields = lines.strip().split('\t')
            offsets[int(fields[0])] = {0:fields[1], 1:fields[2],2:fields[3]}

    # Current support only for quantified read counts from 5' end. Offset from 3' end can be implemented.

    dict_len = {}
    read_count_dict = {}
    # We parse the files for each fragment size which contain the read counts aligned by 5' end for CDS region along with a certain length before and beyond the CDS
    for fsize in range(frag_min, frag_max + 1):
        read_count_dict[fsize] = {}
        with open(infolder+'Read_counts_' + str(fsize) + '.tab') as count_file:
            for lines in count_file:
                fields = lines.strip().split('\t')
                gene = fields[0]
                start_index, length = int(fields[1]), int(fields[2])
                reads_list = map(int, fields[3].split(','))
                dict_len[gene] = length
                read_count_dict[fsize][gene] = {}
                # We start by counting all the read counts listed from the starting index before the start position to similar length beyond stop position
                for i in range(start_index, abs(start_index)+length):
                    # Since there is no zero position, the listed value will shift directly from -1 to 1. So we need to add 1 here to properly index
                    if i >= 0:
                        try:
                            read_count_dict[fsize][gene][i+1] = reads_list[i - start_index]
                        except IndexError:
                            print('INDEX ERROR: length of read list '+str(len(reads_list))+'. Start index: '+str(start_index)+'. Length of gene '+str(length))
                    else:
                        read_count_dict[fsize][gene][i] = reads_list[i - start_index]
    print('Parsed the CDS read counts')

    # Now we generate the A-site profiles according to offsets for specific fragment size and frames
    asite_dict = {}
    for fsize in range(frag_min, frag_max+1):
        for gene in read_count_dict[fsize]:
            if gene not in asite_dict:
                asite_dict[gene] = {}
            asite = []
            for j in range(1, dict_len[gene]+1):
                asite.append(0)
            for pos in sorted(read_count_dict[fsize][gene]):
                # First step is to get the frame of the nucleotide position. Have to careful before the start position
                if pos > 0:
                    frame = pos % 3 - 1  # 0 ,1 ,-1
                if pos < 0:
                    frame = (pos + 1) % 3 - 1  # 0,1,-1
                # the frame 2 is where we get the value of -1
                if frame == -1:
                    frame = 2
                try:
                    offset = int(offsets[fsize][frame])
                # Fsize and frame combinations with ambigious offsets will be of type '15/18' and hence will give value error. They will made offset=0
                except ValueError:
                    offset = 0
                except IndexError:
                    print('IndexError for fsize '+str(fsize)+' and frame '+str(frame)+' and offsets are '+str(offsets))
                    offset = 0

                if offset != 0:
                    # Only those pos before 0 matter when offseted map to a position in CDS
                    if pos < 0 <= pos+offset < len(asite):
                        asite[pos+offset] += read_count_dict[fsize][gene][pos]
                    elif pos > 0 and pos+offset-1 < len(asite):
                        # -1 because the Asite profile is a list with index 0
                        asite[pos+offset-1] += read_count_dict[fsize][gene][pos]
            # This dictionary will store for every gene the A-site profiles for each fragment size
            asite_dict[gene][fsize] = asite

    # Output file for A-site profiles
    asite_file = open('A-site_profiles.tab', 'w')
    asite_table_dict = {}
    for gene in asite_dict:
        # This adds up the read count at every position from all fragment sizes for every gene
        asite_profile = [sum(x) for x in zip(*asite_dict[gene].values())]
        asite_table_dict[gene] = asite_profile
        asite_file.write(gene + '\t' + str(dict_len[gene]) + '\t' + ','.join(map(str, asite_profile)) + '\n')

    # Dumping a pickle object for easy import in downstream analyses
    Pickle.dump(asite_table_dict, open('A-site_profiles.p', 'wb'))
    asite_file.close()



# Using cutoffs for all possible frag sizes
def select_high_cov_genes(folder, frag_min, frag_max, threshold, three_prime, filter_file, include):

    mul_map_file = os.path.join(folder, MULTIPLE_MAPPED_GENE_READ_COUNTS_FILE)
    # Parsing the cds file and getting the length of each gene transcript
    dict_length = {}
    # get the read dict for every gene according to fragment size and at each gene position
    reads_dict = {}
    # A dict containing the total number of reads mapped to a gene
    total_reads = {}
    for fsize in range(frag_min, frag_max + 1):
        if fsize not in reads_dict:
            reads_dict[fsize] = {}
        with open(folder + 'Read_counts_' + str(fsize) + '.tab') as f:
            for lines in f:
                line_list = lines.strip().split('\t')
                gene_name = line_list[0]
                if gene_name not in reads_dict[fsize]:
                    reads_dict[fsize][gene_name] = {}
                start_index = int(line_list[1])
                gene_len = int(line_list[2])
                if gene_name not in dict_length:
                    dict_length[gene_name] = gene_len
                count_list = map(int, line_list[3].split(','))
                if gene_name not in total_reads:
                    total_reads[gene_name] = 0
                idx = start_index
                for val in count_list:
                    if three_prime:
                        # For 3' end alignments, we include mapped reads upto fragment size length beyond the stop codon
                        if 1 <= idx <= gene_len + fsize:
                            reads_dict[fsize][gene_name][idx] = val
                            total_reads[gene_name] += val
                    else:
                        # We will count reads from fsize positions before the start position upto the last nt position of the gene
                        if -fsize <= idx <= gene_len:
                            reads_dict[fsize][gene_name][idx] = val
                            total_reads[gene_name] += val
                    idx += 1
                    if idx == 0:
                        idx += 1

    print('\nParsed read counts for all fragment sizes ')

    log_file = open(output + 'select_genes.log', 'w')
    # Get the number of mul mapped reads to decide whether to delete the gene or not. If a gene has more than 0.1% of reads multiple mapped, we delete it
    mul_map_dict = {}
    mul_map_gene_reads = {}
    with open(mul_map_file) as f:
        for lines in f:
            line_list = lines.strip().split('\t')
            gene_name = line_list[0]
            read_counts = map(int, line_list[1:])
            mul_map_gene_reads[gene_name] = sum(read_counts)
            idx = 0
            for fsize in range(frag_min, frag_max + 1):
                if fsize not in mul_map_dict:
                    mul_map_dict[fsize] = {}
                mul_map_dict[fsize][gene_name] = {0: read_counts[idx], 1: read_counts[idx + 1], 2: read_counts[idx + 2]}
                idx += 3

        mulfile = open(output + 'Genes_multiple_mapped_reads_agg.tab', 'w')
        # List of all genes which have > 1% mul mapped reads and hence will not be considered
        mul_map_genes = []
        for gene in mul_map_gene_reads:
            if mul_map_gene_reads[gene] > 0:
                if gene not in total_reads:
                    log_file.write(gene + ' not in total reads. It must be overlapping gene\n')
                    continue
                else:
                    perc_mul_map = float(mul_map_gene_reads[gene]) * 100 / float(mul_map_gene_reads[gene] + total_reads[gene])
                    mulfile.write(gene + '\t' + str(perc_mul_map) + '%\t' + str(mul_map_gene_reads[gene]) + '\t' + str(total_reads[gene]) + '\n')
                    if perc_mul_map > 1:
                        mul_map_genes.append(gene)
    print('\nParsed the multiple mapped read counts\n')
    dict_gene = {}
    good_genes = {}
    good_genes_mul_map = {}

    # If the IP algorithm has to be applied for a selected subset of genes, they must be loaded from a text file here. If include option is false, these genes will be excluded
    if filter_file:
        filter_genes = True
        select_genes = []
        with open(filter_file) as f:
            for lines in f:
                fields = lines.strip().split('\t')
                gene = fields[0]
                select_genes.append(gene)
        if include:
            print('Using gene list from file ' + filter_file + ' to apply the IP algorithm only for these ' + str(len(select_genes))+' genes')
        else:
            print('Excluding '+str(len(select_genes))+' genes from the application of IP algorithm as provided in the file '+filter_file)
    else:
        filter_genes = False

    print('Starting to select genes for each fragment size and frame')
    for fsize in range(frag_min, frag_max + 1):
        good_genes[fsize] = {}
        good_genes_mul_map[fsize] = {}
        for frame in xrange(3):
            good_genes[fsize][frame] = []
            good_genes_mul_map[fsize][frame] = []
        for gene_name, dict_reads in reads_dict[fsize].iteritems():
            # If we want to exclude or include a set of genes
            if filter_genes:
                # If the include boolean is false, we will exclude the genes in the list select_genes
                if not include and gene_name in select_genes:
                    continue
                # Include only those genes which  are in the list select_genes
                elif include and gene_name not in select_genes:
                    continue
            if gene_name not in dict_gene:
                dict_gene[gene_name] = {}
                for size in range(frag_min, frag_max + 1):
                    dict_gene[gene_name][size] = {}
            reads = []

            try:
                short_utr = False
                # 5' end
                if not three_prime:
                    start = -fsize
                    end = dict_length[gene_name] + 1
                else:
                    start = 1
                    end = dict_length[gene_name] + fsize + 1
                for k in range(start, end):
                    # Ignoring since there is no zero position
                    if k == 0:
                        continue
                    try:
                        reads.append(dict_reads[int(k)])

                    except KeyError:
                        log_file.write(gene_name + '\t' + str(k) + '\tline 257\t' + str(dict_length[gene_name]))
                        short_utr = True
            except KeyError:
                print('Length not available for gene ' + gene_name)
                short_utr = False
            for frame in xrange(3):
                if three_prime:
                    ref = reads[frame::3]

                else:
                    # The extra number of nucleotides on one side of CDS
                    if short_utr:
                        extra_s = min(dict_reads)
                        if extra_s < 12:
                            continue
                    else:
                        extra_s = fsize
                    # and then choose the reads in positions of multiples of 3 according to the frame
                    ref = reads[extra_s % 3 + frame::3]
                    # select the reads of the given 5' end frame  # last_off is -fsize. That many positions will be left out from the right side of the list (after gene ending)
                if gene_name in mul_map_genes:
                    mul_map_ref = mul_map_dict[fsize][gene_name][frame]

                avg_reads = np.mean(ref)

                dict_gene[gene_name][fsize][frame] = (sum(ref), len(ref), avg_reads)
                if avg_reads > threshold:
                    if gene_name in mul_map_genes:
                        try:
                            perc_mul_map_fs_fr = float(mul_map_ref) * 100 / float(mul_map_ref + sum(ref))
                        except ZeroDivisionError:
                            perc_mul_map_fs_fr = 0
                        if perc_mul_map_fs_fr > 1:
                            good_genes_mul_map[fsize][frame].append(gene_name)
                        else:
                            good_genes[fsize][frame].append(gene_name)
                    else:
                        good_genes[fsize][frame].append(gene_name)

    outfile = open(output + 'Sorted_genes_by_Avg_reads_frag_size.tab', 'w')
    for name in dict_gene:
        outfile.write(name)
        for fsize in range(frag_min, frag_max + 1):
            try:
                for frame in xrange(3):
                    outfile.write('\t' + str(dict_gene[name][fsize][frame][0]) + '\t' + str(dict_gene[name][fsize][frame][1]) + '\t' + str(dict_gene[name][fsize][frame][2]))
            except KeyError:
                continue
        outfile.write('\n')

    filtered_cds_dict = {}
    for i in range(frag_min, frag_max + 1):
        filtered_cds_dict[i] = {}
        for frame in xrange(3):
            filtered_cds_dict[i][frame] = {}
    log_file = open(output + 'Gene_filter_statistics.tab', 'w')
    for fsize in filtered_cds_dict:
        for frame in xrange(3):
            for name in good_genes[fsize][frame]:
                filtered_cds_dict[fsize][frame][name] = reads_dict[fsize][name]

            log_file.write(str(fsize) + '\t' + str(frame) + '\t' + str(len(filtered_cds_dict[fsize][frame])) + '\t' + str(len(good_genes_mul_map[fsize][frame])) + '\n')
    return filtered_cds_dict, dict_length


def asite_algorithm_improved_second_offset_correction(reads_dict, 
                                                      dict_len, 
                                                      frag_min, 
                                                      frag_max, 
                                                      output, 
                                                      offset_threshold, 
                                                      off_correction_threshold, 
                                                      three_prime,
                                                      get_coverage_stats=True, 
                                                      advanced=True, 
                                                      conserve_frame=True, 
                                                      bootstrap=False, 
                                                      cov_range=(1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 25, 30, 35, 40, 45, 50)):
    print('Got gene lengths for ' + str(len(dict_len)) + ' genes')

    # Used for debugging
    offset_dic = {}
    correction_dic = {}
    details = {}
    debug_details = {}

    for frame in xrange(3):
        details[frame] = {}
        debug_details[frame] = {}
    log_file = open(output + 'asite_ip.log', 'w')

    # The following dict will contain meta gene for every gene in every fsize and frame. The meta-data include no. of zeros, perc zeros, avg, avg at start and avg at end
    sum_total = {}
    dict_cov_info = {}
    # Initialize files for all fragment size and frame combinations
    for fsize in range(frag_min, frag_max + 1):
        dict_cov_info[fsize] = {}
        sum_total[fsize] = {}
        for frame in xrange(3):
            dict_cov_info[fsize][frame] = {}
            sum_total[fsize][frame] = 0

    # MAIN ANALYSIS STARTS HERE
    for fsize in range(frag_min, frag_max + 1):
        print('Starting analysis for fragment size ' + str(fsize))

        # The following will be used as an index in the read dictionary for each gene
        last_off = -fsize

        for frame in xrange(3):
            if frame not in offset_dic:
                offset_dic[frame] = {}
                correction_dic[frame] = {}
            if fsize not in offset_dic[frame]:
                offset_dic[frame][fsize] = {}
            '''
            FOLLOWING CODE FOR GENES IN EVERY FSIZE AND FRAME
            '''
            for gene in reads_dict[fsize][frame]:
                skip_gene = False
                short_utr = False

                """
                ***  PARSING THE READ VALUES   ***
                """
                # The following will be reads dictionary with nucleotide position as key and number of reads as value
                dict_reads = reads_dict[fsize][frame][gene]

                reads = []
                try:
                    # 5' end
                    if not three_prime:
                        start = -fsize
                        end = dict_len[gene] + 1
                    # For 3' end, we get the reads at the CDS nucleotide positions as well as fsize length after the end of CDS
                    else:
                        start = 1
                        end = dict_len[gene] + fsize + 1

                    for k in range(start, end):
                        # Ignoring since there is no zero position
                        if k == 0:
                            continue
                        try:
                            reads.append(dict_reads[k])
                        except KeyError:
                            log_file.write('ERROR: --Short UTR--KeyError in reads dictionary for ' + gene + ' at position ' + str(k) + 'with gene length of ' +
                                           str(dict_len[gene]) + '\n')
                            short_utr = True
                            if k > 0:
                                skip_gene = True
                                log_file.write('ERROR: KeyError in reads dictionary for ' + gene + ' at position ' + str(k) + 'with gene length of ' + str(dict_len[gene]) + '\n')
                                # Using break instead of continue as this will break the inner for loop and continue the outer for loop
                                break
                except KeyError:
                    # Length not available for this gene
                    skip_gene = True
                    log_file.write('ERROR: KeyError in length dict for ' + gene + '\n')

                if skip_gene:
                    # If gene does not have proper read values, we will not consider it for analysis and hence we remove it and move to next gene in the for loop
                    continue
                # The extra number of nucleotides on one side of CDS
                if short_utr:
                    extra_s = min(dict_reads)
                    if extra_s < 12:
                        continue
                else:
                    extra_s = fsize

                if three_prime:
                    ref = [0] * frame
                    for r in reads[frame::3]:
                        ref += [r, 0, 0]
                else:
                    # To adjust for the len and replace zeroes for out of frame reads, we do the following. This will be followed by deletion/addition of
                    # additional zeroes to make length equal to original seq
                    ref = [0] * (extra_s % 3 + frame)
                    # Choose the reads in positions of multiples of 3 according to the frame
                    for r in reads[extra_s % 3 + frame::3]:  # select the reads of the given 5' end frame
                        ref += [r, 0, 0]
                ref = ref[:-2]  # we exclude the last [0,0] we added at the end
                if (len(reads) - len(ref)) > 0:
                    ref += (len(reads[:last_off]) - len(ref)) * [0]
                # we put it back to the original length (which might have changed when multiple of 3).
                if three_prime:
                    avg_reads = np.mean(ref[frame::3])
                else:
                    avg_reads = np.mean(ref[extra_s % 3 + frame::3])
                """
                ***    CALCULATING THE OFFSET FOR THE GENE  ***
                """
                # try:
                if conserve_frame:
                    if three_prime:
                        best_start_index, best_end_index, offset, score_per_offset = put_reads_in_orf_3_end(ref, extra_s, dict_len[gene], centre_on_p_not_on_a=False,
                                                                                                            advanced=advanced,
                                                                                                            go_three_at_a_time=True)
                    else:
                        best_start_index, best_end_index, offset, score_per_offset = put_reads_in_orf(ref, extra_s, dict_len[gene], centre_on_p_not_on_a=False,
                                                                                                      advanced=advanced, go_three_at_a_time=True)
                else:
                    # we get the offset to centre on frame 1, hence +1 and -2. (-2 is because we use the length to get to the last index by start+len,
                    # but we gave +1 to start so we need to take out an extra one to len).
                    best_start_index, best_end_index, offset, score_per_offset = put_reads_in_orf(ref, extra_s + 1, dict_len[gene] - 2, centre_on_p_not_on_a=False,
                                                                                                  advanced=advanced, go_three_at_a_time=True)

                # If secondary selection criteria is to be applied, we compare the scores of top two offset and the reads near the start codon.
                if advanced:
                    # sort the offsets based on scores.  If scores are same for two or more offsets, they will be sorted according to offset values.
                    sorted_scores = sorted(sorted(score_per_offset), key=score_per_offset.get, reverse=True)
                    # Quality check to make sure the highest offset is the same as the best offset we get from our function
                    if sorted_scores[0] != offset:
                        log_file.write('ERROR: Sorted offsets do not match the offset we get from put_reads_in_orf function for gene ' + gene + ' in fragment size' +
                                       str(fsize) + 'and frame ' + str(frame) + '. The sorted scores are ' + str(sorted_scores) + ' and the offset is ' + str(offset) + '\n')

                    # Difference of top two offsets
                    diff = score_per_offset[offset] - score_per_offset[sorted_scores[1]]

                    # If the difference in scores is less than the avg number of reads across the gene, we apply the secondary selection criteria
                    if diff < avg_reads:
                        # Offsets wit diff less than avg will be listed in here
                        list_offsets_to_compare = []
                        # Add the top two offsets to the list of offsets to compare
                        if score_per_offset[sorted_scores[0]] >= score_per_offset[sorted_scores[1]]:
                            list_offsets_to_compare.append(sorted_scores[0])
                            list_offsets_to_compare.append(sorted_scores[1])
                            # If any other offsets down the order have equal scores with second best offset, then they get added to the list as well
                            for i in range(2, len(sorted_scores)):
                                if score_per_offset[sorted_scores[i]] == score_per_offset[sorted_scores[1]]:
                                    list_offsets_to_compare.append(sorted_scores[i])
                        # The offsets for which the condition is met will be added in here
                        off_true = []
                        # The dict will contain the difference between the average of R2, R3 and R4 and the reads in first codon R1
                        diff_dict = {}

                        # Check the secondary selection criteria of the listed offsets
                        for off in list_offsets_to_compare:
                            # quality check.
                            if off > fsize:
                                log_file.write('Unusual offset ' + str(off) + ' being considered for fsize ' + str(fsize) + ' frame ' + str(frame) + ' in gene ' + gene + '\n')
                                continue
                            # quality check
                            if off % 3 != 0:
                                log_file.write('Unusual offset ' + str(off) + ' being considered for fsize ' + str(fsize) + ' frame ' + str(frame) + ' in gene ' + gene + '\n')
                            # Getting the first 4 codon values in the particular offset
                            if three_prime:
                                reads_offset = reads[off:off + 12]
                            else:
                                reads_offset = reads[extra_s - off:extra_s - off + 12]
                            if not reads_offset:
                                print('Reads offset list is empty')
                                print(extra_s, off)
                            # Checking the condition whether the R1 is less than one-fifth of the average of R2, R3 and R4
                            bool_off, diff_avg = secondary_selection_conditions(reads_offset, frame, threshold=off_correction_threshold)
                            # Adding this offset to the list if the condition is met
                            if bool_off:
                                off_true.append(off)
                                diff_dict[off] = diff_avg
                        # Select the offset which meets the secondary selection criteria
                        if len(off_true) == 1:
                            offset_correction = off_true[0]
                        # If more than one offset meets the secondary selection criteria, then choose the offset with the maximum score
                        elif len(off_true) > 1:
                            diff_compare = {}
                            # check if the scores are equal or not. If the scores are also equal, add the difference of avg(R2, R3, R4) and R1 to diff compare and
                            # choose the offset with max difference
                            max_score = score_per_offset[sorted_scores[0]]
                            for i in range(0, len(off_true)):
                                if score_per_offset[off_true[i]] == max_score:
                                    diff_compare[off_true[i]] = diff_dict[off_true[i]]
                            if len(diff_compare) > 1:
                                sorted_diff = sorted(diff_compare, key=diff_compare.get, reverse=True)
                                offset_correction = sorted_diff[0]
                            else:
                                # check if anything changes
                                # offset_correction = sorted_scores[0]
                                offset_correction = off_true[0]
                        # If no offset meets the secondary selection criteria, we let the offset with the maximum score remain the best offset.
                        # For offsets with equal scores, the smallest offset is the optimal offset
                        else:
                            offset_correction = sorted_scores[0]

                        ''' 
                        CHECK FOR OFFSET CORRECTION 
                        '''
                        # If the offset after applying secondary selection criteria is not the same as offset with top score, change the optimal offset
                        if offset_correction != offset:
                            offset = offset_correction
                        # If the offset is the same as the original one, no need to change the offset
                        elif offset_correction == offset:
                            pass
                        # If the offset correction did not yield any offset (secondary conditions not met for any offset), the initial offset remains the optimal offset
                        elif offset_correction == '':
                            pass

                # OFFSET IS ASSIGNED FOR THE GENE IN THIS FRAGMENT SIZE AND FRAME
                if not skip_gene:
                    offset_dic[frame][fsize][gene] = offset
                else:
                    continue

                if three_prime:
                    sum_total[fsize][frame] += sum(ref[frame::3])
                else:
                    sum_total[fsize][frame] += sum(ref[extra_s % 3 + frame::3])
                dict_cov_info[fsize][frame][gene] = avg_reads

    """
    ***   ANLAYSE DATA FOR DIFFERENT COVERAGE THRESHOLDS ***
    """
    if get_coverage_stats:
        # This dict will contain the bootstrapped distributions from which the error bars will be calculated
        bootstrap_dict = {}
        dict_most_prob_offsets = {}
        print('Running the coverage analysis with offset threshold ' + str(offset_threshold) + '%')
        # Get coverage statistics
        for fsize in range(frag_min, frag_max + 1):
            # Needed to create the error bars for the plot for Figure 3.
            bootstrap_dict[fsize] = {}
            dict_most_prob_offsets[fsize] = {}
            for frame in xrange(3):
                bootstrap_dict[fsize][frame] = {}
                read_avg_dict = {}
                for off in range(0, fsize, 3):
                    bootstrap_dict[fsize][frame][off] = {}
                    read_avg_dict[off] = []
                dict_most_prob_offsets[fsize][frame] = {'off': '', 'perc': ''}
                gene_count = 0
                # Append the meta gene properties to their respective dictionaries according to the offsets
                for gene, offset in offset_dic[frame][fsize].iteritems():
                    try:
                        # Get the meta data of each gene from the dict_cov_info dictionary
                        read_avg = dict_cov_info[fsize][frame][gene]
                    except KeyError:
                        print('KeyError for dict_cov_info in ' + str(fsize) + ' and frame ' + str(frame))
                        continue
                    read_avg_dict[offset].append(read_avg)
                    gene_count += 1
                list_genes = []  # Required for bootstrapping

                for off in range(0, fsize, 3):
                    for j in read_avg_dict[off]:
                        # this is appending the offset value with the read average
                        list_genes.append((off, j))
                if len(list_genes) == 0:
                    continue

                got_offset = False
                # First list will contain the cutoff values and the second list will contain the percentage values
                trend_list = [[], []]
                for c in cov_range:
                    sum_dict, perc_dict = count_stuff_all_offsets(c, read_avg_dict, False)
                    no_of_genes = sum(sum_dict.values())
                    sorted_perc = sorted(perc_dict, key=perc_dict.get, reverse=True)
                    if not got_offset:
                        if no_of_genes >= 10:
                            if perc_dict[sorted_perc[0]] >= offset_threshold:
                                dict_most_prob_offsets[fsize][frame]['off'] = str(sorted_perc[0])
                                dict_most_prob_offsets[fsize][frame]['perc'] = str(perc_dict[sorted_perc[0]])
                                next_highest_perc = perc_dict[sorted_perc[1]]
                                next_best_off = sorted_perc[1]
                                got_offset = True
                            else:
                                dict_most_prob_offsets[fsize][frame]['off'] = str(sorted_perc[0]) + '/' + str(sorted_perc[1])
                                dict_most_prob_offsets[fsize][frame]['perc'] = str(perc_dict[sorted_perc[0]]) + '/' + str(perc_dict[sorted_perc[1]])

                        elif no_of_genes < 10 and dict_most_prob_offsets[fsize][frame]['off'] == '':
                            dict_most_prob_offsets[fsize][frame]['off'] = 'NA'
                            dict_most_prob_offsets[fsize][frame]['perc'] = str(perc_dict[sorted_perc[0]])
                    elif perc_dict[sorted_perc[0]] >= offset_threshold and no_of_genes >= 10:
                        dict_most_prob_offsets[fsize][frame]['perc'] = str(perc_dict[sorted_perc[0]])
                    if no_of_genes >= 10:
                        trend_list[0].append(c)
                        trend_list[1].append(perc_dict[sorted_perc[0]])

                # When the threshold is crossed at an higher cutoff, we do a linear trend analysis to check whether the trend is significant
                if trend_list[1] and trend_list[1][0] < offset_threshold and '/' not in dict_most_prob_offsets[fsize][frame]['off'] and 'NA' not in \
                        dict_most_prob_offsets[fsize][frame]['off']:
                    b1, b0, r, pval, std_err = stat.linregress(trend_list[0], trend_list[1])

                    if pval > 0.05:
                        dict_most_prob_offsets[fsize][frame]['off'] += '/' + str(next_best_off)
                        dict_most_prob_offsets[fsize][frame]['perc'] += '/' + str(next_highest_perc) + '(NS, slope = ' + str(np.round(b1, 3)) + ' pval =' + str(
                            np.round(pval, 3)) + ')'
                    else:
                        dict_most_prob_offsets[fsize][frame]['perc'] += '(Significant, slope = ' + str(np.round(b1, 3)) + ' pval =' + str(np.round(pval, 3)) + ')'

                if bootstrap:
                    mega_sum_dict, mega_perc_dict = bootstrap_gene_count(c, list_genes)
                    Pickle.dump(mega_perc_dict, open('Pickle_dicts/bootstrap_perc_dict.p', 'wb'))
                    Pickle.dump(mega_sum_dict, open('Pickle_dicts/bootstrap_count_dict.p', 'wb'))
                    # Write the properties of the genes at each offset for eack threshold value
                    for off in range(0, fsize, 3):
                        if off not in mega_perc_dict:
                            print(str(off) + ' was not in mega_perc_dict so adding a empty directory for fsize ' + str(fsize) + ' and frame ' + str(frame))
                            mega_perc_dict[off] = [0]
                        try:
                            avg = np.mean(mega_perc_dict[off])
                        except TypeError:
                            avg = 'TE'
                        if c == 1:
                            bootstrap_dict[fsize][frame][off]['avg'] = avg
                    for off in range(0, fsize, 3):
                        try:
                            se_mean = np.std(mega_perc_dict[off])
                        except TypeError:
                            se_mean = 'TE'

                        if c == 1:
                            bootstrap_dict[fsize][frame][off]['se_mean'] = se_mean
                    for off in range(0, fsize, 3):
                        try:
                            low_ci = np.percentile(mega_perc_dict[off], 2.5)
                        except TypeError:
                            low_ci = 'TE'

                        if c == 1:
                            bootstrap_dict[fsize][frame][off]['low_ci'] = low_ci
                    for off in range(0, fsize, 3):
                        try:
                            high_ci = np.percentile(mega_perc_dict[off], 97.5)
                        except TypeError:
                            high_ci = 'TE'
                        if c == 1:
                            bootstrap_dict[fsize][frame][off]['high_ci'] = high_ci
                    for off in range(0, fsize, 3):
                        if off not in mega_sum_dict:
                            print(str(off) + ' was not in mega_sum_dict so adding a empty directory for fsize ' + str(fsize) + ' and frame ' + str(frame))
                            mega_sum_dict[off] = [0]

    """
    ***    WRITE THE RESULTS AND PLOT DISTRIBUTION OF OFFSETS ***
    """
    outfile = open(output + "Results_IP_algorithm.tab", "w")
    perc_file = open(output + "Perc_of_genes_for_all_offsets.tab", "w")

    outfile.write('\n\nMost probable Offsets for Fragment Size and Frame (including coverage data)\n')
    outfile.write('Frag size\tFrame_0\tFrame_1\tFrame_2\n')
    for fsize in range(frag_min, frag_max + 1):
        try:
            outfile.write(str(fsize) + '\t' + str(dict_most_prob_offsets[fsize][0]['off']) + '\t' + str(dict_most_prob_offsets[fsize][1]['off']) + '\t' + str(
                dict_most_prob_offsets[fsize][2]['off']) + '\n')
        except KeyError:
            outfile.write(str(fsize) + '\tNA\tNA\tNA\n')

    outfile.write('\n\nPerc of genes (including coverage data)\n')
    outfile.write('Frag size\tFrame_0\tFrame_1\tFrame_2\n')
    for fsize in range(frag_min, frag_max + 1):
        try:
            outfile.write(str(fsize) + '\t' + str(dict_most_prob_offsets[fsize][0]['perc']) + '\t' + str(dict_most_prob_offsets[fsize][1]['perc']) + '\t' + str(
                dict_most_prob_offsets[fsize][2]['perc']) + '\n')
        except KeyError:
            outfile.write(str(fsize) + '\tNA\tNA\tNA\n')

    outfile.write('\n\n\tNumber of genes\nFrag/Frame\t0\t1\t2\n')
    for fsize in range(frag_min, frag_max + 1):
        outfile.write(str(fsize))
        for frame in xrange(3):
            try:
                outfile.write('\t' + str(len(offset_dic[frame][fsize])))
            except KeyError:
                outfile.write('\tNA')
        outfile.write('\n')

    outfile.write('\n\n\tNumber of reads\nFrag/Frame\t0\t1\t2\n')
    for fsize in range(frag_min, frag_max + 1):
        outfile.write(str(fsize))
        for frame in xrange(3):
            outfile.write('\t' + str(sum_total[fsize][frame]))
        outfile.write('\n')

    perc_file.write('Percentage of genes\nFrag/Frame\t0\t1\t2\n')
    for fsize in range(frag_min, frag_max + 1):
        perc_file.write(str(fsize))
        for frame in xrange(3):
            offset_list = {}
            for off in range(0, fsize, 3):
                offset_list[off] = []
            try:
                for gene, val in offset_dic[frame][fsize].iteritems():
                    if val in offset_list:
                        offset_list[val].append(gene)
                for off in sorted(offset_list):
                    try:
                        perc_file.write('\t' + str(float(len(offset_list[off])) * 100 / float(len(offset_dic[frame][fsize]))))
                    except KeyError:
                        perc_file.write('\tNA')
                    except ZeroDivisionError:
                        perc_file.write('\tNA')
            except KeyError:
                pass
            perc_file.write('\t')
        perc_file.write('\n')
    if bootstrap:
        perc_file.write('Perc of gene bootstrap results (AVERAGE)\nFrag\t0\t3\t6\t9\t12\t15\t18\t21\t0\t3\t6\t9\t12\t15\t18\t21\t0\t3\t6\t9\t12\t15\t18\t21\n')

        for fsize in range(frag_min, frag_max + 1):
            outfile.write(str(fsize))
            for frame in xrange(3):
                try:
                    for off in [0, 3, 6, 9, 12, 15, 18, 21]:
                        perc_file.write('\t' + str(bootstrap_dict[fsize][frame][off]['avg']))
                except KeyError:
                    perc_file.write('\tNA')
                perc_file.write('\t')
            perc_file.write('\n')

        perc_file.write('Perc of gene bootstrap results (SE)\nFrag\t0\t3\t6\t9\t12\t15\t18\t21\t0\t3\t6\t9\t12\t15\t18\t21\t0\t3\t6\t9\t12\t15\t18\t21\n')
        for fsize in range(frag_min, frag_max + 1):
            perc_file.write(str(fsize))
            for frame in xrange(3):
                try:
                    for off in [0, 3, 6, 9, 12, 15, 18, 21]:
                        perc_file.write('\t' + str(bootstrap_dict[fsize][frame][off]['se_mean']))
                except KeyError:
                    perc_file.write('\tNA')
                perc_file.write('\t')
            perc_file.write('\n')
        perc_file.write('Perc of gene bootstrap results (LOW CI)\nFrag\t0\t3\t6\t9\t12\t15\t18\t21\t0\t3\t6\t9\t12\t15\t18\t21\t0\t3\t6\t9\t12\t15\t18\t21\n')
        for fsize in range(frag_min, frag_max + 1):
            perc_file.write(str(fsize))
            for frame in xrange(3):
                try:
                    for off in [0, 3, 6, 9, 12, 15, 18, 21]:
                        perc_file.write('\t' + str(bootstrap_dict[fsize][frame][off]['low_ci']))
                except KeyError:
                    perc_file.write('\tNA')
                perc_file.write('\t')
            perc_file.write('\n')

        perc_file.write('Perc of gene bootstrap results (HIGH CI)\nFrag\t0\t3\t6\t9\t12\t15\t18\t21\t0\t3\t6\t9\t12\t15\t18\t21\t0\t3\t6\t9\t12\t15\t18\t21\n')
        for fsize in range(frag_min, frag_max + 1):
            perc_file.write(str(fsize))
            for frame in xrange(3):
                try:
                    for off in [0, 3, 6, 9, 12, 15, 18, 21]:
                        perc_file.write('\t' + str(bootstrap_dict[fsize][frame][off]['high_ci']))
                except KeyError:
                    perc_file.write('\tNA')
                perc_file.write('\t')
            perc_file.write('\n')

    return dict_most_prob_offsets


def count_stuff_all_offsets(count, dict_count, less=True):
    sum_dict = {}
    perc_dict = {}

    # for offset in [0,3,6,9,12,15,18,21]:

    for offset in dict_count:

        if less:
            sum_dict[offset] = sum(i < count for i in dict_count[offset])
        else:
            sum_dict[offset] = sum(i > count for i in dict_count[offset])

    for offset in dict_count:
        try:
            perc_dict[offset] = np.round(float(sum_dict[offset]) * 100 / float(sum(sum_dict.values())), 2)

        except ZeroDivisionError:
            perc_dict[offset] = 'NA'

    return sum_dict, perc_dict


#  Improved check of secondary offset correction.  Just checks the condition, whether it is valid or not
def secondary_selection_conditions(reads, frame, threshold=5):
    try:
        avg_three_codons = float(reads[3 + frame] + reads[6 + frame] + reads[9 + frame]) / 3
    except IndexError:
        print('IndexError in calculating average of R2,R3,R4')
        print(reads)
        print(frame)
        avg_three_codons = 0
    bool_first = False

    # The average of 2nd, 3rd and 4th codon should be greater than reads in first codon and reads in second codon must be greater in third and fourth codons
    if (avg_three_codons > threshold * reads[frame]) and (reads[3 + frame] > reads[6 + frame]):
        bool_first = True
    diff = avg_three_codons - reads[frame]

    return bool_first, diff


def put_reads_in_orf(reads, gene_start_index, gene_len, centre_on_p_not_on_a=False, advanced=True, go_three_at_a_time=False):
    """
    tries to put most of the reads in the ORF.
    Returns best_start_index, best_end_index, offset (the indices are in reads, offset is how much - in nucleotides - one should move the reads to be aligned with the A site).
     Note that if centre_on_P_not_on_A is False then best_end_index-best_start_index+1 == gene_len-3 since we exclude the first codon when we centre on the A site
     (translation initiate at with AUG in the P site).
     ASSUMES missing reads value is zero
     Assumes the offset is positive (it never moves backward.. reverse genes should be 'straightened beforehand (this happens when they are annotated with reads from parsing of
     SAM file).

    Parameters
    ----------
    go_three_at_a_time
    advanced
    centre_on_p_not_on_a
    gene_len
    gene_start_index
    reads
    """
    if not centre_on_p_not_on_a:
        # Start from second codon
        gene_start_index += 3
        # The length of gene will be counted from second codon until stop codon
        gene_len -= 3
    off = gene_start_index
    # It looks like iteration will happen in opposite direction Starting from the fragment size value+3 (since we are dealing from second codon).
    # However the offset is 0 below even though we are looking at pos with index off and off+gene_len. We are using pos indexes to get sum of the reads while offset is 0.
    reads_in_cds = sum(reads[off: off + gene_len])  # this is for the initial guessed start_index equal to gene_start_index
    score = reads_in_cds
    best_offset = off
    score_per_offset = None
    if advanced:
        score_per_offset = OrderedDict()
        # gene_start_index-off is our actual offset and now we create a dict of scores with offsets as keys
        score_per_offset[gene_start_index - off] = reads_in_cds
    # We are increasing the offset by 3 by shifting the pos indexes by -3 positions to calculate sum of reads
    if go_three_at_a_time:
        off -= 3
    else:
        off -= 1
    while off > 3:
        if go_three_at_a_time:
            reads_in_cds = reads_in_cds + sum(reads[off:off + 3]) - sum(reads[off + gene_len:off + 3 + gene_len])
        else:
            reads_in_cds = reads_in_cds + reads[off] - reads[off + gene_len]
        if advanced:
            score_per_offset[gene_start_index - off] = reads_in_cds
        if reads_in_cds > score:
            # If the read sum is greater than the best score, make it the best score and the offset
            score = reads_in_cds
            best_offset = off
        if go_three_at_a_time:
            off -= 3
        else:
            off -= 1

    return best_offset, best_offset + gene_len - 1, gene_start_index - best_offset, score_per_offset


def put_reads_in_orf_3_end(reads, fsize, gene_len, centre_on_p_not_on_a=False, advanced=True, go_three_at_a_time=False):
    """
    tries to put most of the reads in the ORF.
    Returns best_start_index, best_end_index, offset (the indices are in reads, offset is how much - in nucleotides - one should move the reads to be aligned with the A site).
     Note that if centre_on_P_not_on_A is False then best_end_index-best_start_index+1 == gene_len-3 since we exclude the first codon when we centre on the A site
     (translation initiate at with AUG in the P site).
     ASSUMES missing reads value is zero
     Assumes the offset is positive (it never moves backward.. reverse genes should be 'straightened beforehand (this happens when they are annotated with reads from parsing of
     SAM file).

    Parameters
    ----------
    go_three_at_a_time
    advanced
    centre_on_p_not_on_a
    gene_len
    fsize
    reads
    """
    if not centre_on_p_not_on_a:
        # Start from second codon
        gene_start_index = 3
        # The length of gene will be counted from second codon until stop codon
        gene_len -= 3
    off = gene_start_index
    # Iteration will happen in opposite direction. Starting from the fragment size value+3 (since we are dealing from second codon).
    # However the offset is 0 below even though we are looking at pos with index off and off+gene_len. We are using pos indexes to get sum of the reads while offset is 0.
    reads_in_cds = sum(reads[off: off + gene_len])  # this is for the initial guessed start_index equal to gene_start_index
    score = reads_in_cds
    best_offset = off
    score_per_offset = None
    if advanced:
        score_per_offset = OrderedDict()
        # gene_start_index-off is our actual offset and now we create a dict of scores with offsets as keys
        score_per_offset[gene_start_index - off] = reads_in_cds
    # We are increasing the offset by 3 by shifting the pos indexes by -3 positions to calculate sum of reads
    # if go_three_at_a_time:
    #     off += 3
    # else:
    #     off += 1
    while off < fsize:
        if go_three_at_a_time:
            reads_in_cds = reads_in_cds - sum(reads[off:off + 3]) + sum(reads[off + gene_len:off + 3 + gene_len])
        else:
            reads_in_cds = reads_in_cds + reads[off] - reads[off + gene_len]
        if go_three_at_a_time:
            off += 3
        else:
            off += 1
        if advanced:
            score_per_offset[off - gene_start_index] = reads_in_cds
        if reads_in_cds > score:
            # If the read sum is greater than the best score, make it the best score and the offset
            score = reads_in_cds
            best_offset = off

    return best_offset, best_offset + gene_len - 1, best_offset - gene_start_index, score_per_offset


def bootstrap_gene_count(cutoff, list_genes):
    mega_sum_dict = {}
    mega_perc_dict = {}
    for itr in range(1, 10001):
        # For each iteration, declare a random dict
        rand_dict = {}
        # Select randomly with replacement the indexes of the actual dict
        list_rand_idx = list(np.random.choice(range(0, len(list_genes)), size=len(list_genes)))
        # Populate the random dict with the actual value of offset and the read avg
        for l in list_rand_idx:
            offset, value = list_genes[l]
            if offset not in rand_dict:
                rand_dict[offset] = []
            # Append the gene value
            rand_dict[offset].append(value)
        # Calculate the no of genes for each offset along with the percentage of genes according to offsets for the total genes at this threshold
        sum_dict, perc_dict = count_stuff_all_offsets(cutoff, rand_dict, less=False)
        for off in [0, 3, 6, 9, 12, 15, 18, 21]:
            if off not in mega_sum_dict:
                mega_sum_dict[off] = []
                mega_perc_dict[off] = []
            try:
                mega_sum_dict[off].append(sum_dict[off])
                mega_perc_dict[off].append(perc_dict[off])
            except KeyError:
                mega_sum_dict[off].append(0)
                mega_perc_dict[off].append(0)

    return mega_sum_dict, mega_perc_dict

