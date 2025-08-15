# This is the base dockerfile for the GATK SV pipeline that adds R, a few R packages, and GATK
ARG SAMTOOLS_CLOUD_IMAGE=samtools-cloud:latest
ARG VIRTUAL_ENV_IMAGE=sv-base-virtual-env:latest
ARG GATK_COMMIT="725f9b0ffb1eeac8c9f6da317855cc46b6f137fa"
ARG GATK_JAR="/opt/gatk.jar"
ARG R_INSTALL_PATH=/opt/R


#################### get alias for updating required dependencies
FROM $SAMTOOLS_CLOUD_IMAGE as samtools_cloud


#################### remove unneeded R library files, build GATK, and combine required packages before copying
FROM $VIRTUAL_ENV_IMAGE as virtual_env_image
RUN rm_unneeded_r_library_files.sh

ARG GATK_BUILD_DEP="git git-lfs openjdk-17-jdk"
ARG GATK_RUN_DEP="openjdk-17-jre-headless libgomp1"
ARG GATK_COMMIT
ARG GATK_JAR
ARG DEBIAN_FRONTEND=noninteractive
# install build deps
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-upgrade --no-install-recommends \
                 $GATK_BUILD_DEP $GATK_RUN_DEP
# clone requested GATK_COMMIT
RUN GIT_LFS_SKIP_SMUDGE=1 git clone https://github.com/broadinstitute/gatk.git && \
    cd gatk && \
    git checkout ${GATK_COMMIT}
# build GATK_JAR
RUN cd gatk && \
    ./gradlew localJar && \
    mv $(readlink -f build/libs/gatk.jar) ${GATK_JAR}
# store the required run dependencies, combined from samtools-cloud image, virtual image, and GATK install
RUN mv $APT_REQUIRED_PACKAGES /tmp/r-packages.list
COPY --from=samtools_cloud $APT_REQUIRED_PACKAGES $APT_REQUIRED_PACKAGES
RUN export NEW_PACKAGES=$(diff_of_lists.sh "$(union_lists.sh  /tmp/r-packages.list  "$GATK_RUN_DEP")" \
                                           $APT_REQUIRED_PACKAGES) && \
    printf " $NEW_PACKAGES" | fix_spaces.sh >> $APT_REQUIRED_PACKAGES


################## copy GATK and R installs to final image
FROM $SAMTOOLS_CLOUD_IMAGE

# Copy virtual environment
COPY --from=virtual_env_image /opt /opt
# Set $GATK_JAR
ARG GATK_JAR
ENV GATK_JAR=$GATK_JAR

# Install libraries required for runtime
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-upgrade --no-install-recommends $(fix_spaces.sh $APT_REQUIRED_PACKAGES) && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_*

# Add R executables to PATH.
ARG R_INSTALL_PATH
ARG R_INSTALL_BIN=$R_INSTALL_PATH/bin
ENV PATH=$R_INSTALL_BIN:$PATH

# Test that R and GATK can execute
RUN R --version
RUN java -jar ${GATK_JAR}
