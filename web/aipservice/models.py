from django.db import models
from django.conf import settings

class Job(models.Model):
    INPUT = 'input'
    email = models.EmailField(null=True, blank=True)
    create_date = models.DateTimeField(auto_now_add=True)
    task_id = models.CharField(max_length=256, null=True, blank=True)
    status = models.CharField(max_length=10)
    sam_file = models.FileField(upload_to=INPUT)
    annotation_file = models.FileField(upload_to=INPUT)
    fasta_file = models.FileField(upload_to=INPUT)
    offset_file = models.FileField(null=True, blank=True, upload_to=INPUT)
    filter_file = models.FileField(null=True, blank=True, upload_to=INPUT)
    include = models.BooleanField(null=True)
    min_frag = models.IntegerField()
    max_frag = models.IntegerField()
    three_prime = models.BooleanField()
    overlap = models.IntegerField()
    threshold_avg_reads = models.IntegerField()
    threshold_gene_pct = models.IntegerField()
    threshold_start_codon = models.IntegerField()
    alignment_type = models.CharField(max_length=100)
    get_asite = models.BooleanField()