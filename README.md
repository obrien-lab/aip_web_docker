
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

4. Open a web browser and go to http://0.0.0.0. You can now use the web application.

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

### Reference
[1] Ahmed, N., Sormanni, P., Ciryam, P. et al. Identifying A- and P-site locations on ribosome-protected mRNA fragments using Integer Programming. Sci Rep 9, 6256 (2019). https://doi.org/10.1038/s41598-019-42348-x



