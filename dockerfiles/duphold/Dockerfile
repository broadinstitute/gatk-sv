FROM brentp/musl-hts-nim:latest
RUN apk update && apk upgrade && apk add bash && /bin/bash
RUN cd /  &&      \
    git clone -b master --depth 1 git://github.com/brentp/duphold.git && \
    cd duphold && \
    nimble install -y && \
    nim c -d:release -d:danger -o:/usr/bin/duphold --passC:-flto src/duphold && \
    rm -rf /duphold && duphold --help
ENV GATK_JAR="/opt/gatk.jar"
RUN apk add git git-lfs openjdk8 && \
    cd / && git clone https://github.com/broadinstitute/gatk.git && \
    cd gatk && git lfs install && \
    ./gradlew localJar && \
    mv $(readlink -f build/libs/gatk.jar) ${GATK_JAR} && \
    cd / && \
    java -jar ${GATK_JAR}
