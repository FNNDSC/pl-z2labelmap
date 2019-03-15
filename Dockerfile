# Docker file for the z2labelmap plugin app

FROM fnndsc/ubuntu-python3:latest
MAINTAINER fnndsc "dev@babymri.org"

ENV APPROOT="/usr/src/z2labelmap"  VERSION="0.1"
COPY ["z2labelmap", "${APPROOT}"]
COPY ["requirements.txt", "${APPROOT}"]

WORKDIR $APPROOT

RUN pip install -r requirements.txt

CMD ["z2labelmap.py", "--help"]