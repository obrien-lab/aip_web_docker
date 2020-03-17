from django.urls import path

from . import views

urlpatterns = [
  path('', views.HomeView.as_view(), name='home'),
  path('submit_offset', views.SubmitOffsetView.as_view(), name='submit_offset'),
  path('submit_profile', views.SubmitProfileView.as_view(), name='submit_profile'),
  path('datasets', views.DatasetsView.as_view(), name='datasets'),
  path('offset_report/<str:job_id>/', views.AipReportView.as_view(), name='offset_report'),
  path('profile_report/<str:job_id>/', views.ProfileReportView.as_view(), name='profile_report'),
  path('ajax/get_offset_results/<str:job_id>/', views.get_offset_results, name='get_offset_results'),
  path('ajax/get_job_statistics/', views.get_job_statistics, name='get_job_statistics'),
  path('download/<path:path>/', views.download, name='download')
]