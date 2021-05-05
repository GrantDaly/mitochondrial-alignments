From ubuntu:18.04

RUN apt-get update && apt-get install -y python \
python-pip \
dnsutils

ENV PATH="/opt/sentieon-genomics/bin/:${PATH}"
COPY sentieon-genomics-202010.01 /opt/sentieon-genomics
COPY gen_credentials.py /opt/gen_credentials.py
COPY sentieonStartup.sh /opt/sentieonStartup.sh
COPY alignment-workflow.sh .
COPY requirements.txt .


RUN pip install -r requirements.txt

#ENTRYPOINT ["/bin/bash", "/opt/sentieonStartup.sh;"]
