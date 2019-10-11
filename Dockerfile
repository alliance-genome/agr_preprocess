FROM agrdocker/agr_base_linux_env:latest

WORKDIR /usr/src/app

RUN mkdir output
RUN mkdir download
RUN mkdir download/organism
RUN mkdir download/tmp

ADD requirements.txt .

RUN pip3 install -r requirements.txt

ADD . .

CMD ["python3", "-u", "src/process/combine_interactions_file.py"]
