FROM ubuntu:16.04
RUN apt-get update -y
RUN apt-get install -y python python-pip python-gdal
#install Flask 
RUN pip install -r requirements.txt
СMD python server.py