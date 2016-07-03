FROM ubuntu:16.04
MAINTAINER Andrey Reznik "andrey.reznik.ce@gmail.com"
RUN apt-get update -y
RUN apt-get install -y python python-pip python-gdal
COPY . /app
WORKDIR /app
#install Flask
RUN pip install -r requirements.txt
ENTRYPOINT ["python"]
CMD ["server.py"]
