# This image can run codes from the svtk and sv-pipeline repos
ARG SV_PIPELINE_BASE_IMAGE=sv-pipeline-base:latest
ARG VIRTUAL_ENV_IMAGE=sv-pipeline-virtual-env:latest

# establish alias to copy from later:
FROM $VIRTUAL_ENV_IMAGE as virtual_env_image
RUN rm_unneeded_r_library_files.sh $VIRTUAL_ENV_PATH

FROM $SV_PIPELINE_BASE_IMAGE

# ##################################################
# Copy updated virtual environment from VIRTUAL_ENV_IMAGE
COPY --from=virtual_env_image $VIRTUAL_ENV_PATH $VIRTUAL_ENV_PATH

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
