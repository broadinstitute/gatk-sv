# This is the base dockerfile for the GATK SV pipeline that adds dependencies
# for the sv-pipeline, sv-pipeline-qc, and sv-pipeline-rdtest
ARG SVBASE_IMAGE=sv-base:latest
ARG VIRTUAL_ENV_IMAGE=sv-pipeline-virtual-env:latest

################# get alias for updating dependencies
FROM $SVBASE_IMAGE as base

################ remove unneeded R library files and update required apt dependencies before copying
FROM $VIRTUAL_ENV_IMAGE as virtual_env_image
RUN rm_unneeded_r_library_files.sh

RUN mv $APT_REQUIRED_PACKAGES /tmp/venv-packages.list
COPY --from=base $APT_REQUIRED_PACKAGES $APT_REQUIRED_PACKAGES

RUN export NEW_PACKAGES=$(diff_of_lists.sh /tmp/venv-packages.list $APT_REQUIRED_PACKAGES) && \
    printf " $NEW_PACKAGES" | fix_spaces.sh >> $APT_REQUIRED_PACKAGES

################ start from sv-base
FROM $SVBASE_IMAGE

# Copy in updated virtual environment
COPY --from=virtual_env_image /opt /opt

# NOTE: restore these when we upgrade pysam to latest version. For now: modern HTSLIB prevents pysam from opening
# ENV HTSLIB_LIBRARY_DIR=/opt/samtools/lib
# ENV HTSLIB_INCLUDE_DIR=/opt/samtools/include
# ENV LD_LIBRARY_PATH=$HTSLIB_LIBRARY_DIR:$LD_LIBRARY_PATH

# copy needed resources from git repo
COPY src/sv-pipeline /opt/sv-pipeline
COPY src/WGD /opt/WGD
COPY src/RdTest /opt/RdTest
COPY src/svtk /opt/svtk
COPY src/svtest /opt/svtest
COPY src/svqc /opt/svqc
ENV PATH=/opt/WGD/bin:/opt/bin:$PATH

# install any missing support libraries, along with libraries needed for build
# Install svtk, svtest, svqc
# purge transient (build) deps, leaving only the required packages
ARG BUILD_DEPS="make cmake g++ gcc gfortran"
ARG DEBIAN_FRONTEND=noninteractive
RUN export APT_TRANSIENT_PACKAGES=$(diff_of_lists.sh "$BUILD_DEPS" $APT_REQUIRED_PACKAGES) && \
    apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-upgrade --no-install-recommends $BUILD_DEPS $(fix_spaces.sh $APT_REQUIRED_PACKAGES) && \
    pip install -e /opt/svtk && \
    pip install -e /opt/svtest && \
    pip install -e /opt/svqc && \
    apt-get -qqy remove --purge $APT_TRANSIENT_PACKAGES && \
    apt-get -qqy autoremove --purge && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_*

# test that the installed packages run
RUN svtk -h
RUN svtest -h
RUN svqc -h
RUN plink2 || true
# the command should be installed, but has no help mode, and outputs an error if not used on data
#RUN king

# Compile and test classes from sv-pipeline
# -Compile StitchFragmentedCNVs Java program
# -Compile StitchFragmentedCNVs unit tests
# -Compile VCFParser unit tests
ENV STITCH_JAR="/opt/sv-pipeline/java/build/StitchFragmentedCNVs.jar"
ARG STITCH_UNIT_TEST_JAR="/opt/sv-pipeline/java/build/StitchFragmentedCNVsUnitTest.jar"
ARG VCF_PARSER_UNIT_TEST_JAR="/opt/sv-pipeline/java/build/VCFParserUnitTest.jar"
ARG BUILD_DEPS="openjdk-8-jdk"
ARG DEBIAN_FRONTEND=noninteractive
RUN export APT_TRANSIENT_PACKAGES=$(diff_of_lists.sh "$BUILD_DEPS" $APT_REQUIRED_PACKAGES) && \
    apt-get -qqy update --fix-missing && \
    apt-get -qqy install --no-upgrade --no-install-recommends $BUILD_DEPS && \
    cd /opt/sv-pipeline/java && \
    mkdir -p build/classes && \
    javac -d build/classes org/broadinstitute/svpipeline/StitchFragmentedCNVs.java org/broadinstitute/svpipeline/VCFParser.java && \
    jar cfe build/StitchFragmentedCNVs.jar "org.broadinstitute.svpipeline.StitchFragmentedCNVs" -C build/classes . && \
    rm -rf build/classes/* && \
    javac -d build/classes org/broadinstitute/svpipeline/StitchFragmentedCNVsUnitTest.java org/broadinstitute/svpipeline/StitchFragmentedCNVs.java org/broadinstitute/svpipeline/VCFParser.java && \
    jar cfe build/StitchFragmentedCNVsUnitTest.jar "org.broadinstitute.svpipeline.StitchFragmentedCNVsUnitTest" -C build/classes . && \
    echo "Running StitchFragmentedCNVsUnitTest..." && \
    java -enableassertions -jar $STITCH_UNIT_TEST_JAR && \
    rm -rf build/classes/* $STITCH_UNIT_TEST_JAR && \
    javac -d build/classes org/broadinstitute/svpipeline/VCFParserUnitTest.java org/broadinstitute/svpipeline/VCFParser.java && \
    jar cfe build/VCFParserUnitTest.jar "org.broadinstitute.svpipeline.VCFParserUnitTest" -C build/classes . && \
    echo "Running VCFParserUnitTest..." && \
    java -enableassertions -jar $VCF_PARSER_UNIT_TEST_JAR && \
    rm -rf build/classes/* $VCF_PARSER_UNIT_TEST_JAR && \
    apt-get -qqy remove --purge $APT_TRANSIENT_PACKAGES && \
    apt-get -qqy autoremove --purge && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_*
