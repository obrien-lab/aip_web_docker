FROM python:3.7-slim
RUN apt-get update && apt-get -y upgrade && \
	apt-get install -y build-essential wget \
		libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev libcurl3-dev && \
	apt-get clean && apt-get purge && \
	rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

WORKDIR /usr/src

#Samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 && \
	tar jxf samtools-1.10.tar.bz2 && \
	rm samtools-1.10.tar.bz2 && \
	cd samtools-1.10 && \
	./configure --prefix $(pwd) && \
	make

ENV PATH=${PATH}:/usr/src/samtools-1.10 
RUN python -m pip install --upgrade pip

RUN mkdir -p /opt/app
RUN mkdir -p /opt/app/aip
RUN mkdir -p /files
RUN mkdir -p /files/input

COPY requirements.txt /opt/app/
COPY . /opt/app/aip
WORKDIR /opt/app
RUN python -m pip install -r requirements.txt

