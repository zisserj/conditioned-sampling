FROM movesrwth/storm:stable


RUN apt update && apt-get install -y python3 python3-pip python3.13-venv

ENV VIRTUAL_ENV=/opt/venv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"


RUN pip install numpy scipy matplotlib wheel cython --upgrade setuptools

RUN git clone https://github.com/zisserj/dd.git
RUN git clone https://github.com/zisserj/MCCS.git

ENV DD_FETCH=1 DD_CUDD=1 DD_CUDD_ADD=1
RUN pip install dd/. -vvv --use-pep517 --no-build-isolation

WORKDIR /opt/MCCS
# ENTRYPOINT ["python3", "-c", "import scipy as sc"]


