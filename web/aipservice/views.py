import os
import json
from celery import uuid
from celery import current_app
from django.conf import settings
from django.http import JsonResponse
from django.shortcuts import render, redirect, get_object_or_404
from django.views import View
from django.http import JsonResponse, HttpResponse, Http404
from django.db.models import Count
from django.contrib.sites.models import Site
from django.core.exceptions import PermissionDenied
from .tasks import *
from .forms import *
    
FRAMES = 3

def list_files_in_folder(subfolder):
    dest = os.path.join(settings.MEDIA_ROOT, 'input', subfolder)
    files = []
    if os.path.exists(dest):
        files = os.listdir(dest)
        files.sort()
    return files

class HomeView(View):
    def get(self, request):
        return render(request, 'aipservice/home.html')
    
class UserProfileView(View):
    def get(self, request):
        return render(request, 'aipservice/user_profile.html')    
    
    
class UploadDataView(View):      
    def get(self, request):
        if not request.user.is_authenticated:
            return redirect('account_login')
    
        form = UploadFileForm()
        context = {'form': form, 
                   'default_files': list_files_in_folder('default'),
                   'my_files': list_files_in_folder(os.path.join('users', str(request.user.id))),
                   'max_upload_size': get_max_upload_size(),
                   'max_store_days': settings.MAX_STORE_DAYS
                  }
        return render(request, 'aipservice/datasets.html', context)
                         
    def post(self, request):
        if not request.user.is_authenticated:
            raise PermissionDenied
        form = UploadFileForm(request.POST, request.FILES)
        if form.is_valid():
            file = request.FILES['file']
            dest = os.path.join(settings.MEDIA_ROOT, 'input', 'users', str(request.user.id))
            if not os.path.exists(dest):
                os.mkdir(dest)            
            with open(os.path.join(dest, file.name), 'wb+') as destination:
                for chunk in file.chunks():
                    destination.write(chunk)
            return redirect('datasets')
        else:
            context = {'form': form, 
                       'default_files': list_files_in_folder('default'),
                       'my_files': list_files_in_folder(os.path.join('users', str(request.user.id))),
                       'max_upload_size': get_max_upload_size(),
                       'max_store_days': settings.MAX_STORE_DAYS
                      }
            return render(request, 'aipservice/datasets.html', context)
        
    
class SubmitOffsetView(View):
    def get(self, request):
        if not request.user.is_authenticated:
            return redirect('account_login')
        
        form = AsiteOffsetsJobForm()
        return render(request, 'aipservice/submit_offset.html', { 'form': form })
    
    def post(self, request):
        if not request.user.is_authenticated:
            raise PermissionDenied
        
        form = AsiteOffsetsJobForm(request.POST, request.FILES)
        context = {}        

        if form.is_valid():
            job = form.save(commit=False)
            job.status = "PENDING"
            job.user = request.user
            job.save()
                
            current_site = Site.objects.get_current()
            task = offset_task.delay(current_site.domain,
                                job.id,
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
            
            return redirect('offset_report', job_id = job.id)
        else:
            context['form'] = form
            return render(request, 'aipservice/submit_offset.html', context)
        
class SubmitProfileView(View):
    def get(self, request):
        if not request.user.is_authenticated:
            return redirect('account_login')
        
        form = AsiteProfilesJobForm()
        return render(request, 'aipservice/submit_profile.html', { 'form': form })
    
    def post(self, request):
        if not request.user.is_authenticated:
            raise PermissionDenied
        
        form = AsiteProfilesJobForm(request.POST, request.FILES)
        context = {}        

        if form.is_valid():
            job = form.save(commit=False)
            job.user = request.user
            job.status = "PENDING"
            job.save()
                
            current_site = Site.objects.get_current()
            task = profile_task.delay(current_site.domain,
                                job.id,
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

class OffsetReportView(View):
    def get(self, request, job_id):
        job = get_object_or_404(AsiteOffsetsJob, id=job_id)
        
        if not request.user.is_authenticated:
            return redirect('account_login')
        elif not (request.user.is_superuser or job.user == request.user):
            raise PermissionDenied
            
        folder = os.path.join(settings.MEDIA_ROOT, "Offset_%s" % job_id)
        log_path = os.path.join(folder, "aip.log")
        if not os.path.exists(log_path):
            log_path = None
            
        offset_path = os.path.join(folder, 'A-site_offsets.tab')
        if not os.path.exists(offset_path):
            offset_path = None
        
        perc_gene_path = os.path.join(folder, 'Perc_of_genes_for_all_offsets.tab')
        if not os.path.exists(perc_gene_path):
            perc_gene_path = None

        profile_path = os.path.join(folder, 'A-site_profiles.tab')
        if not os.path.exists(profile_path):
            profile_path = None
        
        return render(request, 'aipservice/offset_report.html', {"job": job, "log_path": log_path, "offset_path": offset_path, "perc_gene_path": perc_gene_path, "profile_path": profile_path})
    
class ProfileReportView(View):
    def get(self, request, job_id):    
        job = get_object_or_404(AsiteProfilesJob, id=job_id)
        
        if not request.user.is_authenticated:
            return redirect('account_login')
        elif not (request.user.is_superuser or job.user == request.user):
            raise PermissionDenied
            
        folder = os.path.join(settings.MEDIA_ROOT, "Profile_%s" % job_id)
        
        log_path = os.path.join(folder, "aip.log")
        if not os.path.exists(log_path):
            log_path = None
        
        profile_path = os.path.join(folder, 'A-site_profiles.tab')
        if not os.path.exists(profile_path):
            profile_path = None
        return render(request, 'aipservice/profile_report.html', {"job": job, "log_path": log_path, "profile_path": profile_path})
    
class JobListView(View):
    def get(self, request):
        if not request.user.is_authenticated:
            return redirect('account_login')
        
        jobs = [{"title": "A-site Offset Jobs", 
                 "job_list": AsiteOffsetsJob.objects.filter(user=request.user),
                 "link": "offset_report"},
                {"title": "A-site Profile Jobs",
                 "job_list": AsiteProfilesJob.objects.filter(user=request.user),
                 "link": "profile_report"}]
        return render(request, 'aipservice/job_list.html', {"jobs": jobs})

def get_job_statistics(request): 
    if not request.user.is_superuser:
        raise PermissionDenied
            
    job_count = Job.objects.count()
    user_count = User.objects.count()
    github_stats = None
    if 'GITHUB_STATS_FILE' in os.environ:
        filepath = os.path.join(settings.MEDIA_ROOT, os.environ['GITHUB_STATS_FILE'])
        if os.path.exists(filepath):
            with open(filepath, 'r') as f:
                github_stats = json.load(f)

    return JsonResponse({'job_count': job_count,
                         'user_count': user_count,
                         'github_stats': github_stats
                        })

def get_offset_results(request, job_id):
    folder  = os.path.join(settings.MEDIA_ROOT, "Offset_%s" % job_id)
    filepath = os.path.join(folder, "A-site_offsets.tab")
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
                    row = line.split('\t')
                    while len(row) < 4:
                        row = row.append("NA")
                    probable_offsets.append(row)
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
    
def delete_file(request, file):
    path = os.path.join(settings.MEDIA_ROOT, 'input', 'users', str(request.user.id), file)
    if os.path.exists(path):
        os.remove(path)
        return redirect('datasets')
    