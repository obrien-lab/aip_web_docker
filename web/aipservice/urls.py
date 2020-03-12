from django.urls import path

from . import views

urlpatterns = [
  path('', views.HomeView.as_view(), name='home'),
  path('submit_aip', views.SubmitAipView.as_view(), name='submit_aip'),
  path('submit_profile', views.SubmitProfileView.as_view(), name='submit_profile'),
  path('datasets', views.DatasetsView.as_view(), name='datasets'),
  path('aip_report/<str:job_id>/', views.AipReportView.as_view(), name='aip_report'),
  path('profile_report/<str:job_id>/', views.ProfileReportView.as_view(), name='profile_report'),
  path('ajax/get_aip_results/<str:job_id>/', views.get_aip_results, name='get_aip_results'),
  path('download/<path:path>/', views.download, name='download')
]