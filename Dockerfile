ARG DOCKER_PULL_TAG=stage
ARG REG=agrdocker
FROM ${REG}/agr_base_linux_env:stage

WORKDIR /usr/src/app

RUN mkdir output
RUN mkdir download
RUN mkdir download/tmp
RUN mkdir download_genetic
RUN mkdir download_genetic/tmp
RUN mkdir download_molecular
RUN mkdir download_molecular/tmp

ADD requirements.txt .

RUN . /root/venv/bin/activate

RUN pip install -r requirements.txt

ADD . .

CMD ["python", "-u", "src/aggregate_preprocessor.py"]
