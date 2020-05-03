from aip import *
# Create your tests here.

folder = "../../../scratch"
bam_file = "../../../aip-files/input/default/Simulated_profiles_23_27_frames_0_1_offset_15.sam"
annotation_file = "../../../aip-files/input/default/Selected_50_genes_23_27_frames_0_1_offset_15.tab" 
fasta_file = "../../../aip-files/input/default/sacCer3_R64-2-1_genome.fa"
offset_file = "../../../aip-files/input/default/A-site_IP_offset_table_yeast.tab"
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
get_profile = True
species = "yeast"

# cdsparser(annotation_file, fasta_file) 

#annotation_file = processAnnotationFile(folder, annotation_file)

run_offset(folder, species, bam_file, annotation_file, fasta_file, min_frag, max_frag, three_prime, overlap, threshold_avg_reads, threshold_gene_pct, threshold_start_codon, filter_file, include, alignment_type, get_profile)


run_profile(folder, species, bam_file, annotation_file, fasta_file, offset_file, min_frag, max_frag, three_prime,  overlap,  alignment_type)