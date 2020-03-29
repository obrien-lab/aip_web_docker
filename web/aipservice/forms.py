import os
from django import forms
from django.conf import settings
from django.template.defaultfilters import filesizeformat
from django.utils.translation import ugettext_lazy as _
from pathlib import Path
from .models import *

def get_max_upload_size():
    root_directory = Path(os.path.join(settings.MEDIA_ROOT, 'input', 'users'))
    current_size = sum(f.stat().st_size for f in root_directory.glob('**/*') if f.is_file() )
    remaining_size = settings.MAX_FOLDER_SIZE - current_size        
    max_upload_size = min(settings.MAX_FILE_SIZE, remaining_size)
    return max_upload_size
        
class UploadFileForm(forms.Form):
    file = forms.FileField()
    
    def clean_file(self):
        file = self.cleaned_data['file']
        max_upload_size = get_max_upload_size()
        
        if max_upload_size <= 0:
            raise forms.ValidationError(_('Not enough space. Please check back later.'))
        elif file.size > max_upload_size:
            raise forms.ValidationError(_('Please keep filesize under %s. Current filesize %s.') % (filesizeformat(max_upload_size), filesizeformat(file.size)))

        return file

class AsiteOffsetsJobForm(forms.ModelForm):
    class Meta:
        model = AsiteOffsetsJob
        fields = ('name', 'species', 'bam_file', 'annotation_file', 'fasta_file', 'filter_file', 'include', 'min_frag', 'max_frag', 'three_prime', 'overlap', 'threshold_avg_reads', 'threshold_gene_pct', 'threshold_start_codon', 'alignment_type', 'get_profile', 'email', )
        
class AsiteProfilesJobForm(forms.ModelForm):
    class Meta:
        model = AsiteProfilesJob
        fields = ('name', 'species', 'bam_file', 'annotation_file', 'fasta_file', 'offset_file', 'min_frag', 'max_frag', 'three_prime', 'overlap', 'alignment_type', 'email', )