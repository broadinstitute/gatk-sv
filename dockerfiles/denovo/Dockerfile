ARG SV_PIPELINE_IMAGE=sv-pipeline:latest
FROM $SV_PIPELINE_IMAGE

COPY /src/denovo /src/denovo/

# make & cpp are needed for installing `tidyverse`.
RUN apt-get -qqy update && \
    apt-get -qqy install ca-certificates build-essential && \
    Rscript -e "\
        install.packages( \
            c( \
                'BiocManager', \
                'Rcpp', \
                'tidyverse', \
                'R.utils', \
                'UpSetR', \
                'ggplotify' \
            ), \
            repos = 'https://cran.rstudio.com', \
            quiet = TRUE \
        )" && \
    pip install fsspec gcsfs
