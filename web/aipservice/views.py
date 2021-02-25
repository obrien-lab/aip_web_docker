import os
import json
import shutil
from datetime import datetime
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
    files = {}
    join = os.path.join
    exists = os.path.exists
    dest = join(settings.MEDIA_ROOT, 'input', subfolder)
    for filetype in ["BAM", "FASTA", "Other"]:
        dest_sub = join(dest, filetype)
        if exists(dest_sub):
            filenames = os.listdir(dest_sub)
            filenames.sort()
            files[filetype] = [{'filename': f, 'filepath': os.path.join(dest_sub, f) } for f in filenames]
    return files


class HomeView(View):
    def get(self, request):
        return render(request, 'aipservice/home.html', {'max_file_size': settings.MAX_FILE_SIZE, 'max_store_days': settings.MAX_STORE_DAYS})
    
    
class UserProfileView(View):
    def get(self, request):
        return render(request, 'aipservice/user_profile.html')    
    
def get_filetype(filename):
    ame, ext = os.path.splitext(filename)
    ext = ext.lower()
    if ext == ".bam" or ext == ".sam":
        filetype = "BAM"
    elif ext == ".fa":
        filetype = "FASTA"
    else:
        filetype = "Other"
    return filetype

class UploadDataView(View):
    def get_context(self, user_id, form):
        context = {'form': form, 
                   'default_files': list_files_in_folder('default'),
                   'max_upload_size': get_max_upload_size(),
                   'max_store_days': settings.MAX_STORE_DAYS
                  }
        
        if user_id:
            my_files = {}
            listfiles = list_files_in_folder(os.path.join('users', user_id))
            for filetype, files in listfiles.items():
                my_files[filetype] = []
                for file in files:
                    my_files[filetype].append({"name": file["filename"], 
                     "path": file["filepath"],
                     "datetime": datetime.datetime.fromtimestamp(os.stat(file["filepath"]).st_mtime),
                     'size': os.path.getsize(file["filepath"])
                    })
            context["my_files"] = my_files
        return context
    
    def get(self, request):    
        form = UploadFileForm()        
        if request.user.is_authenticated:
            user_id = str(request.user.id)
        else:
            user_id = None
        return render(request, 'aipservice/datasets.html', self.get_context(user_id, form))
                         
    def post(self, request):
        if not request.user.is_authenticated:
            raise PermissionDenied
        form = UploadFileForm(request.POST, request.FILES)
        if form.is_valid():
            file = request.FILES['file']
            filetype = get_filetype(file.name)
                
            dest = os.path.join(settings.MEDIA_ROOT, 'input', 'users', str(request.user.id), filetype)
            if not os.path.exists(dest):
                os.makedirs(dest, exist_ok=True)            
            
            with open(os.path.join(dest, file.name), 'wb+') as destination:
                for chunk in file.chunks():
                    destination.write(chunk)
            return redirect('datasets')
        else:
            return render(request, 'aipservice/datasets.html', self.get_context(str(request.user.id), form))
        
    
class SubmitOffsetView(View):
    def get(self, request):        
        form = AsiteOffsetsJobForm()
        return render(request, 'aipservice/submit_offset.html', { 'form': form })
    
    def post(self, request):
        form = AsiteOffsetsJobForm(request.POST, request.FILES)
        context = {}        

        if form.is_valid():
            job = form.save(commit=False)
            job.status = "PENDING"
            if request.user.is_authenticated:
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
            job.task_id = task.id
            job.save()
            return redirect('offset_report', task_id = job.task_id)
        else:
            context['form'] = form
            return render(request, 'aipservice/submit_offset.html', context)
        
class SubmitProfileView(View):
    def get(self, request):        
        form = AsiteProfilesJobForm()
        return render(request, 'aipservice/submit_profile.html', { 'form': form })
    
    def post(self, request):
        form = AsiteProfilesJobForm(request.POST, request.FILES)
        context = {}        

        if form.is_valid():
            job = form.save(commit=False)
            if request.user.is_authenticated:
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
            job.task_id = task.id
            job.save()
            return redirect('profile_report', task_id = job.task_id)
        else:
            context['form'] = form
            return render(request, 'aipservice/submit_profile.html', context)

class OffsetReportView(View):
    def get(self, request, task_id):
        job = get_object_or_404(AsiteOffsetsJob, task_id=task_id)
        species = job.get_species_display()
            
        folder = os.path.join(settings.MEDIA_ROOT, "output", str(job.id))
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
            
        profile_frame0_path = os.path.join(folder, 'A-site_profiles_mapped_to_frame0.tab')
        if not os.path.exists(profile_frame0_path):
            profile_frame0_path = None
        
        return render(request, 'aipservice/offset_report.html', {"job": job, "species": species, "log_path": log_path, "offset_path": offset_path, "perc_gene_path": perc_gene_path, "profile_path": profile_path, "profile_frame0_path": profile_frame0_path})
    
class ProfileReportView(View):
    def get(self, request, task_id):    
        job = get_object_or_404(AsiteProfilesJob, task_id=task_id)
        species = job.get_species_display()
            
        folder = os.path.join(settings.MEDIA_ROOT, "output", str(job.id))
        
        log_path = os.path.join(folder, "aip.log")
        if not os.path.exists(log_path):
            log_path = None
        
        profile_path = os.path.join(folder, 'A-site_profiles.tab')
        if not os.path.exists(profile_path):
            profile_path = None
            
        profile_frame0_path = os.path.join(folder, 'A-site_profiles_mapped_to_frame0.tab')
        if not os.path.exists(profile_frame0_path):
            profile_frame0_path = None
            
        return render(request, 'aipservice/profile_report.html', {"job": job, "species": species, "log_path": log_path, "profile_path": profile_path, "profile_frame0_path": profile_frame0_path})
    
class JobListView(View):
    def get(self, request):
        if not request.user.is_authenticated:
            return redirect('account_login')
        
        jobs = [{"title": "A-site Offset Jobs", 
                 "job_list": AsiteOffsetsJob.objects.filter(user=request.user).order_by('-id'),
                 "link": "offset_report"},
                {"title": "A-site Profile Jobs",
                 "job_list": AsiteProfilesJob.objects.filter(user=request.user).order_by('-id'),
                 "link": "profile_report"}]
        return render(request, 'aipservice/job_list.html', {"jobs": jobs})

def get_job_statistics(request): 
    if not request.user.is_superuser:
        raise PermissionDenied
            
    job_count = Job.objects.count()
    email_count = Job.objects.values("email").distinct().count()
    github_stats = None
    if 'GITHUB_STATS_FILE' in os.environ:
        filepath = os.path.join(settings.MEDIA_ROOT, os.environ['GITHUB_STATS_FILE'])
        if os.path.exists(filepath):
            with open(filepath, 'r') as f:
                github_stats = json.load(f)

    return JsonResponse({'job_count': job_count,
                         'email_count': email_count,
                         'github_stats': github_stats
                        })

def get_offset_results(request, job_id):
    folder = os.path.join(settings.MEDIA_ROOT, "output", job_id)
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
                        row.append("NA")
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
        for line in file:
            if line[0] == "P" or line[0] == "F":
                continue
            blocks = line.strip().split('\t\t')
            for frame in range(FRAMES):
                fields = blocks[frame].split('\t')
                
                if frame == 0:
                    frag = fields[0]
                    fields = fields[1:]
                
                offsets = [x * 3 for x in range(len(fields))]
                offset_pct[frame].append({"x": offsets,
                                          "y": fields,
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
    filetype = get_filetype(file)
    path = os.path.join(settings.MEDIA_ROOT, 'input', 'users', str(request.user.id), filetype, file)
    if os.path.exists(path):
        os.remove(path)
        return redirect('datasets')
    
def cancel_job(request, job_id):
    job = get_object_or_404(Job, id=job_id)

    if not request.user.is_authenticated:
        return redirect('account_login')
    elif not (request.user.is_superuser or job.user == request.user):
        raise PermissionDenied
        
    try:
        current_app.control.revoke(job.task_id, terminate=True)
    except:
        pass

    # clean up the space
    scratch = os.path.join(settings.MEDIA_ROOT, "output", str(job_id), "scratch")
    if os.path.exists(scratch):
         shutil.rmtree(scratch)

    job.status = "CANCELED"
    job.save()
    return redirect('job_list')
    
