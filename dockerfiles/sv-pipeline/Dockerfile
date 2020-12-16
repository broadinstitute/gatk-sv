# This image can run codes from the svtk and sv-pipeline repos

# IMPORTANT: these arguments must be specified at the begining to take advantage of multi-stage build AND runtime specification of base images
ARG SV_PIPELINE_BASE_IMAGE=gatksv/sv-pipeline-base:latest
FROM ${SV_PIPELINE_BASE_IMAGE}

# Install R & bioconductor packages
ARG SV_PIPELINE_BIOCONDUCTOR_PKGS="SNPRelate"
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 make cmake automake && \
    mkdir -p /tmp/R_pkg_download/ && \
    cd /opt/ && \
    Rscript --vanilla install_bioconductor_packages.R ${SV_PIPELINE_BIOCONDUCTOR_PKGS} && \
    cd "/usr/local/lib/R/site-library" && eval ${SLIM_R_LIB_CMD} && \
    cd "/usr/lib/R/site-library" && eval ${SLIM_R_LIB_CMD} && \
    apt-get -qqy remove --purge make cmake automake && \
    apt-get -qqy autoremove --purge && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_* && \
    bcftools --help

# Install plink2 & KING (for relatedness inference)
ARG KING_URL="https://www.kingrelatedness.com/executables/Linux-king215.tar.gz"
ARG PLING2_URL="https://github.com/chrchang/plink-ng/releases/download/2019/plink2_linux_x86_64_20190107.zip"
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 wget unzip && \
    cd /opt && \
    wget -q ${KING_URL} && \
    tar -xzf Linux-king215.tar.gz && \
    wget -q ${PLING2_URL} && \
    unzip plink2_linux_x86_64_20190107.zip && \
    apt-get -qqy purge wget && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_*
ENV PATH="/opt/:${PATH}"

# Compile StitchFragmentedCNVs Java program
ENV STITCH_JAR="/opt/sv-pipeline/java/build/StitchFragmentedCNVs.jar"
RUN cd /opt/sv-pipeline/java && \
    mkdir -p build/classes && \
    javac -d build/classes org/broadinstitute/svpipeline/StitchFragmentedCNVs.java org/broadinstitute/svpipeline/VCFParser.java && \
    jar cfe build/StitchFragmentedCNVs.jar "org.broadinstitute.svpipeline.StitchFragmentedCNVs" -C build/classes . && \
    rm -r build/classes

# Compile StitchFragmentedCNVs unit tests
ENV STITCH_UNIT_TEST_JAR="/opt/sv-pipeline/java/build/StitchFragmentedCNVsUnitTest.jar"
RUN cd /opt/sv-pipeline/java && \
    mkdir -p build/classes && \
    javac -d build/classes org/broadinstitute/svpipeline/StitchFragmentedCNVsUnitTest.java org/broadinstitute/svpipeline/StitchFragmentedCNVs.java org/broadinstitute/svpipeline/VCFParser.java && \
    jar cfe build/StitchFragmentedCNVsUnitTest.jar "org.broadinstitute.svpipeline.StitchFragmentedCNVsUnitTest" -C build/classes . && \
    echo "Running StitchFragmentedCNVsUnitTest..." && \
    java -enableassertions -jar $STITCH_UNIT_TEST_JAR && \
    rm -r build/classes $STITCH_UNIT_TEST_JAR

# Compile VCFParser unit tests
ENV VCF_PARSER_UNIT_TEST_JAR="/opt/sv-pipeline/java/build/VCFParserUnitTest.jar"
RUN cd /opt/sv-pipeline/java && \
    mkdir -p build/classes && \
    javac -d build/classes org/broadinstitute/svpipeline/VCFParserUnitTest.java org/broadinstitute/svpipeline/VCFParser.java && \
    jar cfe build/VCFParserUnitTest.jar "org.broadinstitute.svpipeline.VCFParserUnitTest" -C build/classes . && \
    echo "Running VCFParserUnitTest..." && \
    java -enableassertions -jar $VCF_PARSER_UNIT_TEST_JAR && \
    rm -r build/classes $VCF_PARSER_UNIT_TEST_JAR

# Compile and test CleanVCFPart1 Java program
ENV CLEAN_VCF_PART_1_JAR="/opt/sv-pipeline/java/build/CleanVCFPart1.jar"
RUN cd /opt/sv-pipeline/java && \
    mkdir -p build/classes && \
    javac -d build/classes org/broadinstitute/svpipeline/CleanVCFPart1.java org/broadinstitute/svpipeline/VCFParser.java && \
    jar cfe build/CleanVCFPart1.jar "org.broadinstitute.svpipeline.CleanVCFPart1" -C build/classes . && \
    rm -r build/classes

# Compile and test CleanVCFPart1 unit tests
ENV CLEAN_VCF_PART_1_UNIT_TEST_JAR="/opt/sv-pipeline/java/build/CleanVCFPart1UnitTest.jar"
RUN cd /opt/sv-pipeline/java && \
    mkdir -p build/classes && \
    javac -d build/classes org/broadinstitute/svpipeline/CleanVCFPart1UnitTest.java org/broadinstitute/svpipeline/CleanVCFPart1.java org/broadinstitute/svpipeline/VCFParser.java && \
    jar cfe build/CleanVCFPart1UnitTest.jar "org.broadinstitute.svpipeline.CleanVCFPart1UnitTest" -C build/classes . && \
    echo "Running CleanVCFPart1UnitTest..." && \
    java -enableassertions -jar $CLEAN_VCF_PART_1_UNIT_TEST_JAR && \
    rm -r build/classes $CLEAN_VCF_PART_1_UNIT_TEST_JAR && \
    rm -rf /tmp/* /var/tmp/*
