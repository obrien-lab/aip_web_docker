from django.urls import path

from . import views

urlpatterns = [
  path('', views.HomeView.as_view(), name='home'),
  path('report/<str:task_id>/', views.ReportView.as_view(), name='report'),
  path('ajax/get_results/<str:task_id>/', views.get_results, name='get_results'),
  path('download/<str:path>/', views.download, name='download')
]