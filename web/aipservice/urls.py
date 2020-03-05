from django.urls import path

from . import views

urlpatterns = [
  path('', views.HomeView.as_view(), name='home'),
  path('submit_aip', views.SubmitAipView.as_view(), name='submit_aip'),
  path('submit_profile', views.SubmitProfileView.as_view(), name='submit_profile'),
  path('datasets', views.DatasetsView.as_view(), name='datasets'),
  path('report/<str:job_id>/', views.ReportView.as_view(), name='report'),
  path('ajax/get_results/<str:job_id>/', views.get_results, name='get_results'),
  path('download/<path:path>/', views.download, name='download')
]