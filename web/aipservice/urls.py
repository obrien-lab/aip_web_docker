from django.urls import path

from . import views

urlpatterns = [
  path('', views.HomeView.as_view(), name='home'),
  path('datasets', views.UploadDataView.as_view(), name='datasets'),
  path('submit_offset', views.SubmitOffsetView.as_view(), name='submit_offset'),
  path('submit_profile', views.SubmitProfileView.as_view(), name='submit_profile'),
  path('datasets', views.DatasetsView.as_view(), name='datasets'),
  path('job_list', views.JobListView.as_view(), name='job_list'),
  path('offset_report/<str:job_id>/', views.OffsetReportView.as_view(), name='offset_report'),
  path('profile_report/<str:job_id>/', views.ProfileReportView.as_view(), name='profile_report'),
  path('ajax/get_offset_results/<str:job_id>/', views.get_offset_results, name='get_offset_results'),
  path('ajax/get_job_statistics/', views.get_job_statistics, name='get_job_statistics'),
  path('download/<path:path>/', views.download, name='download'),
  path('delete_file/<path:file>/', views.delete_file, name="delete_file")
]