# Build from the repository root. Reference databases are prepared separately.
# docker build -t neoadjvax-agfusion:1.5.0 -f containers/agfusion/Dockerfile .
FROM python:3.11.11-slim
RUN apt-get update && apt-get install -y --no-install-recommends build-essential pandoc \
    && rm -rf /var/lib/apt/lists/*
RUN python -m pip install --no-cache-dir \
    numpy==1.26.4 pandas==2.2.3 matplotlib==3.9.4 biopython==1.85 \
    pyensembl==2.3.13 gtfparse==2.5.0 mysql-connector-python==9.2.0 agfusion==1.5.0
ENV MPLBACKEND=Agg
RUN agfusion --version && pyensembl --help
