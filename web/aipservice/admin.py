from django.contrib import admin

from .models import *

@admin.register(AsiteOffsetsJob)
class AsiteOffsetsJobAdmin(admin.ModelAdmin):
    list_display = ('pk', 'name', 'task_id', 'user', 'create_date', 'status')
    list_filter = ('user', 'status')

    
@admin.register(AsiteProfilesJob)
class AsiteOffsetsJobAdmin(admin.ModelAdmin):
    list_display = ('pk', 'name', 'task_id', 'user', 'create_date', 'status')
    list_filter = ('user', 'status')