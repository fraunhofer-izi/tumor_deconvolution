FROM python:3.7.4-buster
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH
ENV HOME /code
WORKDIR /code
RUN mkdir -p /deconvs
COPY data/LICENSE.GPL ./LICENSE.GPL
RUN apt-get update --fix-missing && \
    apt-get install -y wget bzip2 ca-certificates \
        libglib2.0-0 libxext6 libsm6 libxrender1 \
        git mercurial subversion libopenblas-dev \
        liblapack-dev gfortran inotify-tools rsync
ENV CONDA https://repo.anaconda.com/miniconda/Miniconda3-4.7.12.1-Linux-x86_64.sh
RUN wget --quiet "$CONDA" -O miniconda.sh && \
    bash miniconda.sh -b -p /opt/conda
RUN ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh
RUN . /opt/conda/etc/profile.d/conda.sh && conda activate base
RUN conda update conda
RUN conda install -c conda-forge pymc3
RUN pip install git+https://github.com/pymc-devs/pymc3
COPY data/header ./header
COPY data/grains ./grains
COPY input_test /input_test
COPY data/chars.pkl ./chars.pkl
ENV FEATURES HUGO
COPY data/main.sh ./main.sh
RUN chmod a+x main.sh
COPY data/deconvolute.py ./deconvolute.py
COPY data/stickbreaking.py ./stickbreaking.py
COPY data/warmup.sh ./warmup.sh
RUN ./warmup.sh
ENV MODE normal
ENTRYPOINT ["./main.sh"]
