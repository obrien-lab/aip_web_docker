
### Introduction

This is a containerized web application that runs the Integer Programming method to identify the A-site positions within ribosome protected fragments (A-site IP) [1]. It uses docker and docker-compose for smooth portablity across platforms, including Linux, macOS and Windows. 

The Docker Compose creates and runs five services and each service is in its own container: (1) a web application based on the Python web framework Django; (2) message broker Redis; (3) Python task queuing package Celery; (4) database PostgresSQL; (5) web server Nginx.

### Prerequisite

Make sure that Docker and Docker Compose are installed. You can simply install the desktop version for [Mac](https://docs.docker.com/docker-for-mac/install/) and [Windows](https://docs.docker.com/docker-for-windows/install/). And Docker Compose is included as part of those desktop installs. 

On Linux systems, you need to first install [Docker](https://docs.docker.com/install/) and then [Docker Compose](https://docs.docker.com/compose/install/) separately.

### Quick Start

1. Create a working directory

```
$ mkdir workspace
$ cd workspace

```

2. Clone the code

```
$ git clone https://github.com/obrien-lab/aip_web_docker.git
$ cd aip_web_docker

```

3. Run

```
$ docker-compose up

```

4. Open a web browser and go to http://0.0.0.0. A superuser has already been created with both username and password being "admin". You can now use the web application.


### Configurations

1. Change compose file. By default, "docker-compose up" will use docker_compose.yml as the compose file. You can change to another compose file by adding the "-f" flag. For example. "docker-compose -f production.yml up" will use production.yml as the compose file. 

2. Change [Volume](https://docs.docker.com/storage/volumes/). Docker can manage the files through Volume, and those files can live beyond the lifecycle of the Docker container. The volumes are configured in the compose file. By default, the input and output files of the A-site IP algorithm are mounted to the folder aip-file, which is at the same level of aip_web_docker, e.g.

```
    volumes:
      - ../aip-files:/files

```
You can change "../aip-files" to your preferred location. 


3. Change environment variables. The file .env stores environment variables, such as postgres username and password. You can change it for better security.

4. Change Nginx configurations. You can configure the Nginx server using the file nginx/sites-enabled/django_project. 

5. Add email support. The application can send an email notification to the user once a job has finished. To enable this function, put your email host server information at the bottom of file web/aip_project/settings.py. You also need to add your email account information in the file .env, e.g.

```
EMAIL_HOST_USER=myusername
EMAIL_HOST_PASSWORD=mypassword

```

6. Change the default superuser. Go to the file web/start-server.sh or web/start-server-production.sh to make the changes.

7. Configue authentication providers. This application uses [django-allauth](https://django-allauth.readthedocs.io/en/latest/) for user authentication. Extra configuration may be needed on the provider side. For example, to use Google account authentication, you need to register the application with Google first. Please follow the [intructions](https://django-allauth.readthedocs.io/en/latest/providers.html).

8. Enable password for Redis. By default, Redis password is disabled as in file docker_compose.yml. You can enable Redis password by adding "command: redis-server --requirepass ${REDIS_PASS}" to the Redis service as in file production.yml, and then define the environment variable REDIS_PASS in file .env.

### Advanced usage
The key piece of code to identify A-site offset and generate A-site density profiles is in script https://github.com/obrien-lab/aip_web_docker/blob/master/web/aipservice/aip.py. Specifically, the function "run_offset" is the entry point to identiry A-site offsets and the function"run_profile" is the entry point to generate A-site density profiles. Computationally literate users can run this script in command line or adapt the script to their own need. 

### Reference
[1] Ahmed, N., Sormanni, P., Ciryam, P. et al. Identifying A- and P-site locations on ribosome-protected mRNA fragments using Integer Programming. Sci Rep 9, 6256 (2019). https://doi.org/10.1038/s41598-019-42348-x



