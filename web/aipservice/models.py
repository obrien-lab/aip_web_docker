from django.db import models
from django.conf import settings

class Job(models.Model):
    email = models.EmailField(null=True, blank=True)
    create_date = models.DateTimeField(auto_now_add=True)
    finish_date = models.DateTimeField(null=True)
    task_id = models.CharField(max_length=256, null=True, blank=True)
    status = models.CharField(max_length=10)
    
class AipJob(Job):  
    INPUT = '/files/input'
    SPECIES_LIST = (
        ('yeast', 'S. cerevisiae'),
        ('ecoli', 'E. coli'),
        ('mouse', 'Mouse'),
        ('other', 'Other')
    )
    species = models.CharField(max_length=10, choices=SPECIES_LIST)
    bam_file = models.FilePathField(path = INPUT)
    annotation_file = models.FilePathField(path = INPUT)
    fasta_file = models.FilePathField(path = INPUT)
    filter_file = models.FilePathField(path = INPUT, null=True, blank=True)
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
    
class ProfileJob(Job):    
    INPUT = '/files/input'
    SPECIES_LIST = (
        ('yeast', 'S. cerevisiae'),
        ('ecoli', 'E. coli'),
        ('mouse', 'Mouse'),
        ('other', 'Other')
    )
    species = models.CharField(max_length=10, choices=SPECIES_LIST)
    bam_file = models.FilePathField(path = INPUT)
    annotation_file = models.FilePathField(path = INPUT)
    fasta_file = models.FilePathField(path = INPUT)
    offset_file = models.FilePathField(path = INPUT)
    min_frag = models.IntegerField()
    max_frag = models.IntegerField()
    three_prime = models.BooleanField()
    overlap = models.IntegerField()
    alignment_type = models.CharField(max_length=100)
