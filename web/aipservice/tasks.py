import os
from django.conf import settings
from celery import shared_task
from .models import Job

@shared_task
def aip_task(job_id,
            bam_file, 
            annotation_file, 
            fasta_file,
            offset_file,
            min_frag = 20, 
            max_frag = 35, 
            three_prime = False, 
            overlap = 0, 
            threshold_avg_reads = 1,
            threshold_gene_pct = 70,
            threshold_start_codon = 5,
            filter_file = '',
            include = True,
            alignment_type = "genome",
            get_asite = False
            ):  
    job = Job.objects.get(id = job_id)
    job.status = "RUNNING"
    job.task_id = aip_task.request.id
    job.save()

    # First convert BAM file to SAM format. Requires samtools to be installed.
    folder = os.path.join(settings.MEDIA_ROOT, job_id)
    if not os.path.exists(folder):
        os.makedirs(folder)
    sam_file = processBamFile(folder, bam_file)
    
    '''
    annotation_file = processAnnotationFile(folder, annotation_file)
    
    if alignment_type == "genome" :
        count_dict, mul_count_dict = samparser_genome(sam_file, min_frag, max_frag, three_prime)
        print('Parsed the SAM file. Starting to quantify CDS read counts')
        create_cds_counts_genome(annotation_file, fasta_file, folder, count_dict, mul_count_dict, min_frag, max_frag, three_prime, overlap)
    else:
        count_dict, mul_count_dict, total_dict = samparser_transcriptome(sam_file, min_frag, max_frag, three_prime)
        print('Parsed the SAM file aligned to the transcriptome. Starting to quantify CDS read counts')
        create_cds_counts_transcriptome(annotation_file, fasta_file, folder, count_dict, mul_count_dict, total_dict, min_frag, max_frag, three_prime)

    if get_asite:
        generate_asite_profiles(min_frag, max_frag, offset_file, folder)

    # Filter genes which have greater than threshold (default=1) reads per codon on average
    filtered_genes, dataset_gene_len = select_high_cov_genes(folder, min_frag, max_frag, threshold_avg_reads, three_prime, filter_file, include)
    print('Parsed the CDS file')

    offset_dict = asite_algorithm_improved_second_offset_correction(filtered_genes, dataset_gene_len, min_frag, max_frag, folder, threshold_start_codon, threshold_gene_pct, three_prime)
    '''
    job = Job.objects.get(id = job_id)
    job.status = "SUCCESS"
    job.save()
