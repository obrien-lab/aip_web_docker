from django.db import models
from django.conf import settings
from django.contrib.auth.models import User

class Job(models.Model):
    name = models.CharField(max_length=256)
    user = models.ForeignKey(User, on_delete=models.CASCADE, null=True)
    email = models.EmailField()
    create_date = models.DateTimeField(auto_now_add=True)
    finish_date = models.DateTimeField(null=True)
    task_id = models.CharField(max_length=256, null=True, blank=True)
    status = models.CharField(max_length=10)
    
class AsiteOffsetsJob(Job):  
    SPECIES_LIST = (
        ('yeast', 'S. cerevisiae'),
        ('ecoli', 'E. coli'),
        ('mouse', 'Mouse'),
        ('other', 'Other')
    )
    species = models.CharField(max_length=10, choices=SPECIES_LIST)
    bam_file = models.CharField(max_length=256)
    annotation_file = models.CharField(max_length=256)
    fasta_file = models.CharField(max_length=256)
    filter_file = models.CharField(max_length=256, null=True, blank=True)
    include = models.BooleanField()
    min_frag = models.IntegerField()
    max_frag = models.IntegerField()
    three_prime = models.BooleanField()
    overlap = models.IntegerField()
    threshold_avg_reads = models.IntegerField()
    threshold_gene_pct = models.IntegerField()
    threshold_start_codon = models.IntegerField()
    alignment_type = models.CharField(max_length=100)
    get_profile = models.BooleanField()
    
class AsiteProfilesJob(Job):    
    SPECIES_LIST = (
        ('yeast', 'S. cerevisiae'),
        ('ecoli', 'E. coli'),
        ('mouse', 'Mouse'),
        ('other', 'Other')
    )
    species = models.CharField(max_length=10, choices=SPECIES_LIST)
    bam_file = models.CharField(max_length=256)
    annotation_file = models.CharField(max_length=256)
    fasta_file = models.CharField(max_length=256)
    offset_file = models.CharField(max_length=256)
    min_frag = models.IntegerField()
    max_frag = models.IntegerField()
    three_prime = models.BooleanField()
    overlap = models.IntegerField()
    alignment_type = models.CharField(max_length=100)
