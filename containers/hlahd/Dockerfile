# HLA-HD (Kawaguchi et al. 2017, Human Mutation) for HLA_TYPING.
#
# HLA-HD is free for academic use but distributed only after registration,
# so the tarball is not downloaded here: put your copy next to this file
# and pass its name as HLAHD_TARBALL. See README.md.
#
#   docker build --build-arg HLAHD_TARBALL=hlahd.1.7.1.tar.gz -t neoadjvax/hlahd:1.7.1 containers/hlahd

FROM ubuntu:22.04

ARG HLAHD_TARBALL

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update \
    && apt-get install -y --no-install-recommends bowtie2 g++ make python3 procps \
    && rm -rf /var/lib/apt/lists/*

COPY ${HLAHD_TARBALL} /tmp/hlahd.tar.gz
RUN mkdir -p /opt/hlahd \
    && tar -xzf /tmp/hlahd.tar.gz -C /opt/hlahd --strip-components=1 \
    && rm /tmp/hlahd.tar.gz \
    && cd /opt/hlahd \
    && sh install.sh \
    && test -x /opt/hlahd/bin/hlahd.sh \
    && test -d /opt/hlahd/dictionary \
    && test -d /opt/hlahd/freq_data \
    && test -f /opt/hlahd/HLA_gene.split.txt

ENV PATH=/opt/hlahd/bin:${PATH}
