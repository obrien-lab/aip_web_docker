import os
from celery import uuid
from celery import current_app
from django import forms
from django.conf import settings
from django.http import JsonResponse
from django.shortcuts import render, redirect, get_object_or_404
from django.views import View
from django.http import JsonResponse, HttpResponse, Http404
from .tasks import *
from .models import *
    
FRAMES = 3

class AipJobForm(forms.ModelForm):
    class Meta:
        model = AipJob
        fields = ('species', 'bam_file', 'annotation_file', 'fasta_file', 'filter_file', 'include', 'min_frag', 'max_frag', 'three_prime', 'overlap', 'threshold_avg_reads', 'threshold_gene_pct', 'threshold_start_codon', 'alignment_type', 'get_profile', 'email', )
        
class ProfileJobForm(forms.ModelForm):
    class Meta:
        model = ProfileJob
        fields = ('species', 'bam_file', 'annotation_file', 'fasta_file', 'min_frag', 'max_frag', 'three_prime', 'overlap', 'alignment_type', 'email', )

class HomeView(View):
    def get(self, request):
        return render(request, 'aipservice/home.html')
        
class SubmitAipView(View):
    def get(self, request):
        form = AipJobForm()
        return render(request, 'aipservice/submit_aip.html', { 'form': form })
    
    def post(self, request):
        form = AipJobForm(request.POST, request.FILES)
        context = {}        

        if form.is_valid():
            job = form.save(commit=False)
            job.status = "PENDING"
            job.save()
                
            task = aip_task.delay(job.id,
                                job.species,
                                job.bam_file, 
                                job.annotation_file, 
                                job.fasta_file,
                                job.min_frag, 
                                job.max_frag, 
                                job.three_prime, 
                                job.overlap, 
                                job.threshold_avg_reads,
                                job.threshold_gene_pct,
                                job.threshold_start_codon,
                                job.filter_file,
                                job.include,
                                job.alignment_type,
                                job.get_profile)
            
            return redirect('aip_report', job_id = job.id)
        else:
            context['form'] = form
            return render(request, 'aipservice/submit_aip.html', context)
        
class SubmitProfileView(View):
    def get(self, request):
        form = ProfileJobForm()
        return render(request, 'aipservice/submit_profile.html', { 'form': form })
    
    def post(self, request):
        form = ProfileJobForm(request.POST, request.FILES)
        context = {}        

        if form.is_valid():
            job = form.save(commit=False)
            job.status = "PENDING"
            job.save()
                
            task = profile_task.delay(job.id,
                                job.species,
                                job.bam_file, 
                                job.annotation_file, 
                                job.fasta_file,
                                job.offset_file,
                                job.min_frag, 
                                job.max_frag, 
                                job.three_prime, 
                                job.overlap, 
                                job.alignment_type)
            
            return redirect('profile_report', job_id = job.id)
        else:
            context['form'] = form
            return render(request, 'aipservice/submit_profile.html', context)

class DatasetsView(View):
    def get(self, request):
        return render(request, 'aipservice/datasets.html')

class AipReportView(View):
    def get(self, request, job_id):
        job = get_object_or_404(AipJob, id=job_id)
        folder = os.path.join(settings.MEDIA_ROOT, "AIP_%s" % job_id)
        log_path = os.path.join(folder, "aip.log")
        if not os.path.exists(log_path):
            log_path = None
            
        offset_path = os.path.join(folder, 'A-site_offsets.tab')
        if not os.path.exists(offset_path):
            offset_path = None
        
        profile_path = os.path.join(folder, 'A-site_profiles.tab')
        if not os.path.exists(profile_path):
            profile_path = None
        
        return render(request, 'aipservice/aip_report.html', {"job": job, "log_path": log_path, "offset_path": offset_path, "profile_path": profile_path})
    
class ProfileReportView(View):
    def get(self, request, job_id):
        job = get_object_or_404(ProfileJob, id=job_id)
        folder = os.path.join(settings.MEDIA_ROOT, "Profile_%s" % job_id)
        
        log_path = os.path.join(folder, "aip.log")
        if not os.path.exists(log_path):
            log_path = None
        
        profile_path = os.path.join(folder, 'A-site_profiles.tab')
        if not os.path.exists(profile_path):
            profile_path = None
        return render(request, 'aipservice/profile_report.html', {"job": job, "log_path": log_path, "profile_path": profile_path})
    
def get_aip_results(request, job_id):
    folder  = os.path.join(settings.MEDIA_ROOT, job_id)
    filepath = os.path.join(folder, "Results_IP_algorithm.tab")
    block = ""
    
    # offsets for each fragment size and frame 
    OFFSETS = "Most probable Offsets for Fragment Size and Frame (including coverage data)"
    probable_offsets = []
    
    # number of genes for each fragment size and frame 
    NUMBER_OF_GENES = "Number of genes"
    genes = []
    
    # distribution of mRNA fragment size
    NUMBER_OF_READS = "Number of reads"
    frag_size_hist = []    
    for frame in range(FRAMES):
        frag_size_hist.append({"x": [],
                               "y": [],
                               "type": 'bar',
                               "name": "Frame " + str(frame)
        });
        
    total_reads = 0
    
    # reading the file
    with open(filepath, 'r') as file:
        for line in file:
            line = line.strip()
            
            # offsets
            if line == OFFSETS:
                block = OFFSETS
            elif block == OFFSETS:
                if len(line) == 0:
                    block = ""
                elif line[0] != "F":
                    probable_offsets.append(line.split('\t'))
            else: 
                # genes distribution
                if line == NUMBER_OF_GENES:
                    block = NUMBER_OF_GENES
                    continue
                elif block == NUMBER_OF_GENES: 
                    if len(line) == 0:
                        block = ""
                    elif line[0] != "F":
                        genes.append(line.split('\t'))
                else: 
                    # reads distribution
                    if line == NUMBER_OF_READS:
                        block = NUMBER_OF_READS
                        continue
                    elif block == NUMBER_OF_READS:
                        if len(line) == 0:
                            block = ""
                        elif line[0] != "F":
                            fields = line.split('\t')
                            for frame in range(FRAMES):
                                frag_size_hist[frame]["x"].append(fields[0])
                                reads = int(fields[frame + 1])
                                frag_size_hist[frame]["y"].append(reads)
                                total_reads = total_reads + reads
    for frame in range(FRAMES):
        frag_size_hist[frame]["y"] = [reads/total_reads*100 for reads in frag_size_hist[frame]["y"] ]
    
    # distribution of offset value
    filepath = os.path.join(folder, "Perc_of_genes_for_all_offsets.tab")
    offset_pct = {}
    for frame in range(FRAMES):
        offset_pct[frame] = []
    with open(filepath, 'r') as file:
        for idx, line in enumerate(file):
            fields = line.strip().split('\t')
            frag = fields[0]
            length = int((len(fields) - 1) / FRAMES)
            offsets = list(range(0, length * FRAMES + 1, FRAMES))
            for frame in range(FRAMES):
                offset_pct[frame].append({"x": offsets,
                                          "y": fields[1 + length * frame : 1 + length * frame + length],
                                          "type": 'bar',
                                          "name": frag
                                         }) 
                
    return JsonResponse({'offsets': probable_offsets,
                         'frag_size_hist': frag_size_hist,
                         'offset_pct': offset_pct,
                        'genes': genes
                        })

def download(request, path):
    if os.path.exists(path):
        with open(path, 'rb') as fh:
            response = HttpResponse(fh.read(), content_type="application/vnd.ms-excel")
            response['Content-Disposition'] = 'inline; filename=' + os.path.basename(path)
            return response
    raise Http404