import os
import logging
import datetime
from django.conf import settings
from django.core.mail import send_mail
from celery import shared_task
from .models import *
from .aip import *

logger = logging.getLogger(__name__)

@shared_task
def offset_task(job_id,
            species,
            bam_file, 
            annotation_file, 
            fasta_file,
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
    job = AsiteOffsetsJob.objects.get(id = job_id)
    job.status = "RUNNING"
    job.task_id = offset_task.request.id
    job.finish_date = datetime.datetime.now()
    job.save()

    # create the working folder
    folder = os.path.join(settings.MEDIA_ROOT, "Offset_%d" % job_id)
    if not os.path.exists(folder):
        os.makedirs(folder)

    status = "SUCCESS"
    try:
        run_offset(folder, 
                species,
                bam_file, 
                annotation_file, 
                fasta_file,
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
    job = AsiteOffsetsJob.objects.get(id = job_id)
    job.status = status
    job.save()
    
    # send notification email
    subject = 'A-site IP job finished'
    message = 'Your A-site IP job is finished.'
    email_from = settings.EMAIL_HOST_USER
    recipient_list = [job.email,]
    send_mail( subject, message, email_from, recipient_list )

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
    job = AsiteProfilesJob.objects.get(id = job_id)
    job.status = "RUNNING"
    job.task_id = profile_task.request.id
    job.finish_date = datetime.datetime.now()
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
    job = AsiteProfilesJob.objects.get(id = job_id)
    job.status = status
    job.save()
    