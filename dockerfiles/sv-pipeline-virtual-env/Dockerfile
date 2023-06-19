ARG SV_BASE_MINI_IMAGE=:sv-base-mini:latest
ARG R_VIRTUAL_ENV_IMAGE=sv-base-virtual-env:latest
ARG PYTHON_VIRTUAL_ENV_IMAGE=samtools-cloud-virtual-env:latest

ARG CONDA_INSTALL_DIR=/opt/conda
ARG CONDA_ENV_NAME="gatk-sv"

############# get handle for sv-base-mini, for copying HTSLIB info
FROM $SV_BASE_MINI_IMAGE as sv_base_mini

############# get python environment and add to it
FROM $PYTHON_VIRTUAL_ENV_IMAGE as python_virtual_env

# install any build dependencies
ARG BUILD_DEPS="g++ make apt-transport-https ca-certificates wget libcurl4-openssl-dev"
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-install-recommends $BUILD_DEPS

# install conda packages
# NOTE: need to use scipy=1.7.3 instead of scipy=1.8.0 because it makes hail angry
ARG CONDA_PKGS="cython=0.29.28 numpy=1.22.3 pandas=1.4.2 scipy=1.7.3 scikit-learn=1.0.2 intervaltree=3.1.0 \
                matplotlib=3.5.1 natsort=8.1.0 google-cloud-dataproc=4.0.2 seaborn=0.12.2"
RUN mamba install -qy --freeze-installed -n $CONDA_ENV_NAME -c conda-forge -c bioconda $CONDA_PKGS

# copy in HTSLIB install so that pysam uses same version as is available in pipeline
# NOTE: restore these when we update pysam version. For now, we want to use the older version of HTSLIB that comes with
# pysam, because the newer versions are what introduce the error
#COPY --from=sv_base_mini /opt/samtools /opt/samtools
#ENV HTSLIB_LIBRARY_DIR=/opt/samtools/lib
#ENV HTSLIB_INCLUDE_DIR=/opt/samtools/include
#ENV LD_LIBRARY_PATH=$HTSLIB_LIBRARY_DIR:$LD_LIBRARY_PATH

# As of March 29, 2022 the newest pysam version is 0.19.0, however that version does not allow setting END < POS, which
# causes a problem for our BND handling. The newest version that is compatible with our BND processing is 0.15.4
# In general, pysam needs to be pip installed or built, because you want to use the installed HTSLIB rather than copy in
# duplicate libraries. The older pysam that needs to be installed for BND issues has problems with cython during pip
# install, so must be built from source. One last pysam wrinkle: it's install can't handle the newest versions of
# setuptools, so we're going to downgrade then re-upgrade them just for pysam
ARG PYSAM_VERSION=0.15.4
RUN export SETUPTOOLS_VERSION=$(python -c 'import setuptools; print(setuptools.__version__)') && \
    pip install setuptools==57.5.0 && \
    wget -q https://github.com/pysam-developers/pysam/archive/refs/tags/v$PYSAM_VERSION.tar.gz && \
    tar -xzf v$PYSAM_VERSION.tar.gz && \
    cd pysam-$PYSAM_VERSION && \
    python setup.py build && \
    python setup.py install && \
    pip install setuptools==$SETUPTOOLS_VERSION

# pybedtools needs to be installed via pip because it doesn't like the updated python
# hail's latest version is only available via pip or local build
ARG PIP_PKGS="pybedtools==0.9.0 hail==0.2.93"
RUN pip3 --no-cache-dir install $PIP_PKGS

# clean unneeded stuff
RUN conda clean -ay --force-pkgs-dirs

############## merge envs by copying python virtual env (smaller env) into R virtual env (bigger env)
FROM $R_VIRTUAL_ENV_IMAGE as build_virtual_env

# copy in python env (not the whole conda install: we don't need conda anymore) and set ENV variables
ARG CONDA_INSTALL_DIR
ARG CONDA_ENV_NAME
ARG CONDA_ENV_PATH=$CONDA_INSTALL_DIR/envs/$CONDA_ENV_NAME
COPY --from=python_virtual_env $CONDA_ENV_PATH $CONDA_ENV_PATH
ARG CONDA_ENV_BIN=$CONDA_ENV_PATH/bin
ENV PATH=$CONDA_ENV_BIN:$PATH
ENV CONDA_ENV_NAME=$CONDA_ENV_NAME

# update the lists of transient (needed for install) and support (needed for runtime) packages
ARG BUILD_DEPS="make cmake g++ gcc gfortran \
                libbz2-dev libblas-dev libicu-dev liblapack-dev liblzma-dev libpcre2-dev zlib1g-dev \
                pkg-config libcurl4-openssl-dev libssl-dev libatlas-base-dev libssh2-1-dev libssh2-1 \
                libgit2-dev libfftw3-dev libjpeg8-dev libxml2-dev libfontconfig1-dev libharfbuzz-dev \
                libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev"
ARG RUN_DEPS="file bzip2 libssl3 libatlas3-base libssh2-1 libgit2-1.1 libfftw3-3 libjpeg8 libxml2"
RUN export NEW_PACKAGES=$(diff_of_lists.sh "$RUN_DEPS" $APT_REQUIRED_PACKAGES) && \
    printf " $NEW_PACKAGES" | fix_spaces.sh >> $APT_REQUIRED_PACKAGES

# install R packages
ARG R_PACKAGES="assertthat beeswarm BH BSDA caret cli crayon DAAG data.table devtools digest dplyr e1071 fansi fpc \
                generics gert glue HardyWeinberg hash latticeExtra magrittr Matrix metap mnormt nlme nloptr nnet \
                numDeriv perm pillar pkgconfig plogr plyr purrr pwr R6 RColorBrewer Rcpp reshape reshape2 rlang ROCR \
                rpart stringi stringr survival tibble tidyr tidyselect utf8 vioplot withr zoo"
ARG BIOCONDUCTOR_PKGS="SNPRelate multtest"
ARG DEBIAN_FRONTEND=noninteractive
RUN export APT_TRANSIENT_PACKAGES=$(diff_of_lists.sh "$BUILD_DEPS" $APT_REQUIRED_PACKAGES) && \
    apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-install-recommends $BUILD_DEPS $(fix_spaces.sh $APT_REQUIRED_PACKAGES) && \
    install_bioconductor_packages.R $BIOCONDUCTOR_PKGS && \
    install_deprecated_R_package.sh "https://cran.r-project.org/src/contrib/Archive/MASS/MASS_7.3-58.tar.gz" && \
    install_R_packages.R $R_PACKAGES && \
    apt-get -qqy remove --purge $APT_TRANSIENT_PACKAGES && \
    apt-get -qqy autoremove --purge && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_*

# Install plink2 & KING (for relatedness inference)
ARG KING_URL="https://www.kingrelatedness.com/executables/Linux-king215.tar.gz"
RUN mkdir -p /opt/bin && \
    cd /opt/bin && \
    wget -q $KING_URL && \
    tar -xzf Linux-king215.tar.gz && \
    rm -f Linux-king215.tar.gz

ARG PLINK2_URL="https://github.com/chrchang/plink-ng/releases/download/2019/plink2_linux_x86_64_20190107.zip"
RUN cd /opt/bin && \
    wget -q $PLINK2_URL && \
    unzip plink2_linux_x86_64_20190107.zip && \
    rm -f plink2_linux_x86_64_20190107.zip

ENV PATH=/opt/bin:$PATH
