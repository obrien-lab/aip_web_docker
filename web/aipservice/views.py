import os

from celery import uuid
from celery import current_app
from django import forms
from django.conf import settings
from django.http import JsonResponse
from django.shortcuts import render, redirect, get_object_or_404
from django.views import View
from django.http import JsonResponse, HttpResponse, Http404
from .tasks import aip_task
from .models import Job
    
FRAMES = 3

class JobForm(forms.ModelForm):
    class Meta:
        model = Job
        fields = ('sam_file', 'annotation_file', 'fasta_file', 'offset_file', 'filter_file', 'include', 'min_frag', 'max_frag', 'three_prime', 'overlap', 'threshold_avg_reads', 'threshold_gene_pct', 'threshold_start_codon', 'alignment_type', 'get_asite', 'email', )
    
class HomeView(View):
    def get(self, request):
        form = JobForm()
        return render(request, 'aipservice/home.html', { 'form': form })
    
    def post(self, request):
        form = JobForm(request.POST, request.FILES)
        context = {}        

        if form.is_valid():
            job = form.save(commit=False)
            job.status = "PENDING"
            job.save()
            
            offset_file = job.offset_file.path if job.offset_file else None
            
            filter_file = job.filter_file.path if job.filter_file else None 
                
            task = aip_task.delay(job.id,
                                job.sam_file.path, 
                                job.annotation_file.path, 
                                job.fasta_file.path,
                                offset_file,
                                job.min_frag, 
                                job.max_frag, 
                                job.three_prime, 
                                job.overlap, 
                                job.threshold_avg_reads,
                                job.threshold_gene_pct,
                                job.threshold_start_codon,
                                filter_file,
                                job.include,
                                job.alignment_type,
                                job.get_asite)
            
            return redirect('report', job_id = job.id)
        else:
            context['form'] = form
            return render(request, 'aipservice/home.html', context)
    
class ReportView(View):
    def get(self, request, job_id):
        job = get_object_or_404(Job, id=job_id)
        return render(request, 'aipservice/report.html', {"job": job})
    
def get_results(request, task_id):
    folder  = settings.MEDIA_ROOT
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
    file_path = os.path.join(settings.MEDIA_ROOT, path)
    if os.path.exists(file_path):
        with open(file_path, 'rb') as fh:
            response = HttpResponse(fh.read(), content_type="application/vnd.ms-excel")
            response['Content-Disposition'] = 'inline; filename=' + os.path.basename(file_path)
            return response
    raise Http404