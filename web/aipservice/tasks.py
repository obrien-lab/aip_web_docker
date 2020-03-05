import os
from django.conf import settings
from celery import shared_task
from .models import Job
from .aip import *

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

    # create the working folder
    folder = os.path.join(settings.MEDIA_ROOT, str(job_id))
    if not os.path.exists(folder):
        os.makedirs(folder)

    status = "SUCCESS"
    try:
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
    except Exception as e:
        print("Error running the Asite-IP job: ", e)
        status = "ERROR"
    
    # update the job status
    job = Job.objects.get(id = job_id)
    job.status = status
    job.save()
    

@shared_task
def profile_task(job_id,
            bam_file, 
            annotation_file, 
            fasta_file,
            offset_file,
            min_frag = 20, 
            max_frag = 35, 
            three_prime = False, 
            overlap = 0, 
            alignment_type = "genome",
            ):  
    job = Job.objects.get(id = job_id)
    job.status = "RUNNING"
    job.task_id = profile_task.request.id
    job.save()

    # create the working folder
    folder = os.path.join(settings.MEDIA_ROOT, str(job_id))
    if not os.path.exists(folder):
        os.makedirs(folder)

    status = "SUCCESS"
    try:
        run_profile(folder, 
                bam_file, 
                annotation_file, 
                fasta_file,
                offset_file,
                min_frag, 
                max_frag, 
                three_prime, 
                overlap, 
                alignment_type)
    except Exception as e:
        print("Error getting A-site profiles: ", e)
        status = "ERROR"
    
    # update the job status
    job = Job.objects.get(id = job_id)
    job.status = status
    job.save()
    