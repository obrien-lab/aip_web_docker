import os
import logging
import datetime
from django.conf import settings
from django.core.mail import send_mail
from celery import shared_task
from billiard.exceptions import Terminated
from .models import *
from .aip import *

logger = logging.getLogger(__name__)

def send_notification_mail(job, domain, job_type):
    # send notification email
    if settings.EMAIL_HOST_USER:
        subject = 'A-site IP job finished'
        message = 'Your A-site IP job has finished. Please visit %s/%s_report/%s to view the results.' % (domain, job_type, job.task_id)
        email_from = settings.EMAIL_HOST_USER
        recipient_list = [job.email,]
        send_mail( subject, message, email_from, recipient_list )
    
@shared_task(soft_time_limit=18000, time_limit=18060, throws=(Terminated,))
def offset_task(domain,
            job_id,
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
    try:
        job = AsiteOffsetsJob.objects.get(id = job_id)
        
        if job.status == "CANCELED":
            return
        
        job.status = "RUNNING"
        job.task_id = offset_task.request.id
        job.save()

        # create the working folder
        folder = os.path.join(settings.MEDIA_ROOT, "output", str(job_id))
        if not os.path.exists(folder):
            os.makedirs(folder, exist_ok=True)

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
        status = "SUCCESS"
    except Exception as e:
        logger.error("Error getting Asite-IP offset: %s", str(e))
        status = "ERROR"
    
    # update the job status
    job = AsiteOffsetsJob.objects.get(id = job_id)
    job.status = status
    job.finish_date = datetime.datetime.now()
    job.save()
    send_notification_mail(job, domain, "offset")

@shared_task(soft_time_limit=18000, time_limit=18060, throws=(Terminated,))
def profile_task(domain,
            job_id,
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
    try:
        job = AsiteProfilesJob.objects.get(id = job_id)
        
        if job.status == "CANCELED":
            return
    
        job.status = "RUNNING"
        job.task_id = profile_task.request.id
        job.save()

        # create the working folder
        folder = os.path.join(settings.MEDIA_ROOT, "output", str(job_id))
        if not os.path.exists(folder):
            os.makedirs(folder, exist_ok=True)

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
        status = "SUCCESS"
    except Exception as e:
        logger.error("Error getting A-site profiles: %s", str(e))
        status = "ERROR"
    
    # update the job status
    job = AsiteProfilesJob.objects.get(id = job_id)
    job.status = status
    job.finish_date = datetime.datetime.now()
    job.save()
    send_notification_mail(job, domain, "profile")