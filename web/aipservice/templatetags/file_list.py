import os
from django import template
from django.conf import settings

register = template.Library()

def list_files_in_folder(subfolder):
    dest = os.path.join(settings.MEDIA_ROOT, 'input', subfolder)
    files = []
    if os.path.exists(dest):
        files = os.listdir(dest)
    return [{'filename': file, 'filepath': os.path.join(dest, file) } for file in files]

@register.inclusion_tag('file_list.html')
def file_list(user):    
    return {'default_files': list_files_in_folder('default'), 
            'my_files': list_files_in_folder(os.path.join('users', str(user.id)))}