ARG R_VERSION=4.4.1
FROM r-base:${R_VERSION}


RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        curl \
        python3 \
        bedtools \
        python-is-python3 \
        gpg \
    && rm -rf /var/lib/apt/lists/*

RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list \
    && curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | gpg --dearmor -o /usr/share/keyrings/cloud.google.gpg \
    && apt-get update \
    && apt-get install -y --no-install-recommends google-cloud-sdk \
    && rm -rf /var/lib/apt/lists/*

ENTRYPOINT [ "bash" ]
CMD [ "-c" ]

COPY src/calculate_APS.R /src/
COPY src/correct_permutated_SVs.R /src/
COPY src/generate_SV_permutations.R /src/
COPY src/generate_cwas_metrics.sh /src/
COPY src/integrate_SV_vs_gencode.R /src/
COPY src/integrate_stats_across_nc_elements.R /src/
COPY src/reorganize_sv_vs_conservative.R /src/

