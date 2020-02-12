#!/usr/bin/env bash
# start-server.sh
cd /opt/app/aip

python manage.py migrate 

echo "from django.contrib.auth.models import User; User.objects.count()>0 or User.objects.create_superuser('admin', 'admin@example.com', 'admin')" | python manage.py shell

/usr/local/bin/gunicorn aip_project.wsgi:application -w 2 -b :8000
