FROM python:2.7

RUN apt-get update && apt-get -y install apt-file && apt-file update && apt-get -y install vim

RUN mkdir /queue && mkdir /queue/public && mkdir /queue/templates && mkdir /app

COPY public /app/public
COPY templates /app/templates
COPY server-queue.py /app
COPY requirements.txt /app
COPY filters.txt /app

WORKDIR /app

RUN pip install -r /app/requirements.txt

#CMD /bin/bash
CMD /usr/local/bin/python server-queue.py