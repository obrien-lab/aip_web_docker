from aip import *
# Create your tests here.

folder = "../../../scratch"
bam_file = "../../../aip-files/sacCer3_Pop.bam"
annotation_file = "../../../aip-files/sacCer3_CDS.tab" 
fasta_file = "../../../aip-files/sacCer3_R64-2-1_genome.fa"
offset_file = ""
min_frag = 20 
max_frag = 35 
three_prime = False 
overlap = 0 
threshold_avg_reads = 1
threshold_gene_pct = 70
threshold_start_codon = 5
filter_file = ""
include = True
alignment_type = "genome"
get_asite = False


run_aip(folder, 
        bam_file, 
        annotation_file, 
        fasta_file,
        offset_file,
        min_frag, 
        max_frag, 
        three_prime, 
        overlap, 
        threshold_avg_reads,
        threshold_gene_pct,
        threshold_start_codon,
        filter_file,
        include,
        alignment_type,
        get_asite) 
