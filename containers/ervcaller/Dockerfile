# Build on a Linux Docker host; transfer the image to Gadi as a SIF.
# ERVcaller is CC BY-NC 4.0. TE reference libraries have separate terms.
FROM ubuntu:22.04
ARG DEBIAN_FRONTEND=noninteractive
ARG ERVCALLER_COMMIT=802dbf49ef61e722f16f84c0a6feb1283d8a360a
LABEL org.opencontainers.image.source="https://github.com/xunchen85/ERVcaller"
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates curl build-essential zlib1g-dev libbz2-dev liblzma-dev \
    bwa samtools r-base-core perl python3 \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /opt/ERVcaller
RUN curl -fsSL "https://github.com/xunchen85/ERVcaller/archive/${ERVCALLER_COMMIT}.tar.gz" \
      | tar xz --strip-components=1 \
    && tar xzf Scripts/SE-MEI-master.tar.gz -C /opt \
    && make -C /opt/SE-MEI extractSoftclipped \
    && install /opt/SE-MEI/extractSoftclipped /usr/local/bin/extractSoftclipped \
    && perl -c /opt/ERVcaller/ERVcaller_v1.4.pl \
    && printf '%s\n' "${ERVCALLER_COMMIT}" > /opt/ERVcaller/UPSTREAM_COMMIT
ENV PATH="/opt/ERVcaller:/usr/local/bin:${PATH}"
WORKDIR /work
