import os
import logging
from django.conf import settings
from celery import shared_task
from .models import *
from .aip import *

logger = logging.getLogger(__name__)

@shared_task
def aip_task(job_id,
            species,
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
            get_profile = False
            ):  
    job = AipJob.objects.get(id = job_id)
    job.status = "RUNNING"
    job.task_id = aip_task.request.id
    job.save()

    # create the working folder
    folder = os.path.join(settings.MEDIA_ROOT, "AIP_%d" % job_id)
    if not os.path.exists(folder):
        os.makedirs(folder)

    status = "SUCCESS"
    try:
        run_aip(folder, 
                species,
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
                get_profile)
    except Exception as e:
        logger.error("Error getting Asite-IP offset: %s", str(e))
        status = "ERROR"
    
    # update the job status
    job = AipJob.objects.get(id = job_id)
    job.status = status
    job.save()
    

@shared_task
def profile_task(job_id,
            species,
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
    job = ProfileJob.objects.get(id = job_id)
    job.status = "RUNNING"
    job.task_id = profile_task.request.id
    job.save()

    # create the working folder
    folder = os.path.join(settings.MEDIA_ROOT, "Profile_%d" % job_id)
    if not os.path.exists(folder):
        os.makedirs(folder)

    status = "SUCCESS"
    try:
        run_profile(folder, 
                species,
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
        logger.error("Error getting A-site profiles: %s", str(e))
        status = "ERROR"
    
    # update the job status
    job = ProfileJob.objects.get(id = job_id)
    job.status = status
    job.save()
    