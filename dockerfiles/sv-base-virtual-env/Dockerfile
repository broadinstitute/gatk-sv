ARG UBUNTU_RELEASE="18.04"
ARG R_RELEASE_VERSION="3.6.3"
ARG GATK_COMMIT="a33bf19dd3188af0af1bd17bce015eb20ba73227"
ARG GATK_JAR="/opt/gatk.jar"
ARG VIRTUAL_ENV_PATH="/gatk-sv-env"

ARG R_INSTALL_BIN=/opt/R/$R_RELEASE_VERSION/bin

####################################################################################################

# install GATK (needed for PrintSVEvidence)
FROM ubuntu:$UBUNTU_RELEASE as build_gatk
ARG GATK_TRANSIENT_DEP="git git-lfs"
ARG GATK_DEP="openjdk-8-jdk ${GATK_TRANSIENT_DEP}"
ARG GATK_COMMIT
ARG GATK_JAR
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-upgrade --no-install-recommends \
                 ${GATK_DEP} && \
    cd / && git clone https://github.com/broadinstitute/gatk.git && \
    cd gatk && git lfs install && \
    git checkout ${GATK_COMMIT} && ./gradlew localJar && \
    mv $(readlink -f build/libs/gatk.jar) ${GATK_JAR}

##################################################

FROM ubuntu:$UBUNTU_RELEASE as build_virtual_env
ARG SV_BASE_R_PKGS="optparse BiocManager"
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8
ARG CONDA_INSTALL_DIR
ARG CONDA_RELEASE
ARG VIRTUAL_ENV_PATH
ARG DEBIAN_FRONTEND=noninteractive

# get prerequisites installing R
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 make git wget curl apt-transport-https ca-certificates gnupg-agent software-properties-common gdebi-core

# install R
ARG R_RELEASE_VERSION
RUN curl -O https://cdn.rstudio.com/r/ubuntu-1804/pkgs/r-${R_RELEASE_VERSION}_1_amd64.deb && \
    gdebi -n -q r-${R_RELEASE_VERSION}_1_amd64.deb
# get the list of packages needed to support R
RUN dpkg -I r-3.6.3_1_amd64.deb | grep Depends | sed -e 's/Depends://' -e 's/,//g' -e 's/^ *//' \
    >  /opt/r_support_packages.list

ARG R_INSTALL_BIN
ENV PATH=/opt/R/$R_RELEASE_VERSION/bin:$PATH

RUN echo $PATH && ls /opt/R/$R_RELEASE_VERSION/bin/

COPY dockerfiles/sv-base-virtual-env/install_R_packages.R $R_INSTALL_BIN
COPY dockerfiles/sv-base-virtual-env/install_bioconductor_packages.R $R_INSTALL_BIN
COPY dockerfiles/sv-base-virtual-env/install_deprecated_R_package.sh $R_INSTALL_BIN
COPY dockerfiles/sv-base-virtual-env/rm_unneeded_r_library_files.sh $R_INSTALL_BIN

RUN Rscript -e "update.packages(contriburl=contrib.url(repos=c('http://lib.stat.cmu.edu/R/CRAN/', 'https://cran.rstudio.com')), clean=TRUE, quiet=TRUE, ask=FALSE)"

RUN install_R_packages.R $SV_BASE_R_PKGS

####################################################################################################

# Copy python packages from build_virtual_env. Effectively removes everything but the minimal environment and flattens
# docker layers in the final image.
# This virtual environment is "active" automatically as long as $VIRTUAL_ENV_PATH and $VIRTUAL_ENV_PATH/bin is in $PATH
FROM ubuntu:$UBUNTU_RELEASE
ARG VIRTUAL_ENV_PATH
ARG GATK_JAR
ARG R_INSTALL_BIN
ARG R_RELEASE_VERSION
ARG DEBIAN_FRONTEND=noninteractive

# make the virtual environment path
RUN echo $VIRTUAL_ENV_PATH && mkdir -p $VIRTUAL_ENV_PATH/bin

COPY --from=build_virtual_env /opt/r_support_packages.list $VIRTUAL_ENV_PATH/

# ensure this image is up-to-date
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
          make git wget curl apt-transport-https ca-certificates gnupg-agent software-properties-common && \
    xargs apt-get -qqy install --no-install-recommends < $VIRTUAL_ENV_PATH/r_support_packages.list && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_*

# copy GATK jar
COPY --from=build_gatk $GATK_JAR $VIRTUAL_ENV_PATH
ENV GATK_JAR=$VIRTUAL_ENV_PATH/gatk.jar

# copy R install
RUN mkdir -p $VIRTUAL_ENV_PATH/R/
COPY --from=build_virtual_env /opt/R/$R_RELEASE_VERSION $VIRTUAL_ENV_PATH/R/$R_RELEASE_VERSION
RUN mkdir -p /opt/R/ &&  ln -s $VIRTUAL_ENV_PATH/R/$R_RELEASE_VERSION /opt/R/$R_RELEASE_VERSION

ENV R_RELEASE_VERSION=$R_RELEASE_VERSION
ENV VIRTUAL_ENV_PATH=$VIRTUAL_ENV_PATH
ENV PATH=$VIRTUAL_ENV_PATH:$VIRTUAL_ENV_PATH/bin:/opt/R/$R_RELEASE_VERSION/bin:$PATH
