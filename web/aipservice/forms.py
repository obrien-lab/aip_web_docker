from django import forms
from .models import *

class UploadFileForm(forms.Form):
    file = forms.FileField()

class AsiteOffsetsJobForm(forms.ModelForm):
    class Meta:
        model = AsiteOffsetsJob
        fields = ('species', 'bam_file', 'annotation_file', 'fasta_file', 'filter_file', 'include', 'min_frag', 'max_frag', 'three_prime', 'overlap', 'threshold_avg_reads', 'threshold_gene_pct', 'threshold_start_codon', 'alignment_type', 'get_profile', 'email', )
        
class AsiteProfilesJobForm(forms.ModelForm):
    class Meta:
        model = AsiteProfilesJob
        fields = ('species', 'bam_file', 'annotation_file', 'fasta_file', 'offset_file', 'min_frag', 'max_frag', 'three_prime', 'overlap', 'alignment_type', 'email', )