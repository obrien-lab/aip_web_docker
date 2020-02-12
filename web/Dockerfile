FROM python:3.7-slim

RUN python -m pip install --upgrade pip

RUN mkdir -p /opt/app
RUN mkdir -p /opt/app/aip

COPY requirements.txt /opt/app/
COPY . /opt/app/aip
WORKDIR /opt/app
RUN python -m pip install -r requirements.txt
