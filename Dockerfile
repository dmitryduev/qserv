FROM python:2.7

RUN apt-get update && apt-get -y install apt-file && apt-file update && \
    apt-get -y install vim git wget sudo gcc curl

RUN mkdir /queue && mkdir /app

# Create a non-root user roboao and switch to it
RUN adduser --disabled-password --gecos '' --shell /bin/bash roboao && \
    chown -R roboao:roboao /app
RUN echo "roboao ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers.d/90-roboao
USER roboao

RUN sudo chown -R roboao:roboao /queue /app

# All users can use /home/roboao as their home directory
ENV HOME=/home/roboao
RUN chmod 777 /home/roboao

# Install latest Miniconda
RUN curl -so ~/miniconda.sh https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh && \
    chmod +x ~/miniconda.sh && \
    ~/miniconda.sh -b -p ~/miniconda
ENV PATH=/home/roboao/miniconda/bin:$PATH
ENV CONDA_AUTO_UPDATE_CONDA=false

COPY public /app/public
COPY templates /app/templates
COPY server-queue.py /app
COPY requirements.txt /app
COPY filters.txt /app

WORKDIR /app

RUN pip install -r /app/requirements.txt

#CMD /bin/bash
CMD python server-queue.py