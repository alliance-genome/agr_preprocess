ARG ALLIANCE_RELEASE=latest
ARG REG=agrdocker
FROM ${REG}/agr_base_linux_env:${ALLIANCE_RELEASE}

WORKDIR /usr/src/app

RUN mkdir output
RUN mkdir download
RUN mkdir download/tmp
RUN mkdir download_genetic
RUN mkdir download_genetic/tmp
RUN mkdir download_molecular
RUN mkdir download_molecular/tmp

ADD requirements.txt .

RUN pip3 install -r requirements.txt

ADD . .

CMD ["python3", "-u", "src/aggregate_preprocessor.py"]
