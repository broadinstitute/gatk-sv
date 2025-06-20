FROM australia-southeast1-docker.pkg.dev/cpg-common/images/sv/manta:2023-09-14-v0.28.3-beta-3f22f94d

ARG VERSION=1.6.0
ARG MANTA_INSTALL_DIR=/usr/local/bin/manta/
RUN sed -i 's/enableRemoteReadRetrievalForInsertionsInGermlineCallingModes = 1/enableRemoteReadRetrievalForInsertionsInGermlineCallingModes = 0/g' ${MANTA_INSTALL_DIR}/bin/configManta.py.ini
