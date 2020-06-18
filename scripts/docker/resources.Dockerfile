# This docker file serves a single purpose: to store the dependencies in a common place 
# so that other Dockerfiles can fetch needed dependencies.
# The image built from this docker file is intermediary and never meant to be pushed

FROM alpine:latest

ENV RESOURCES="/GATKSV_PIPELINE_V1_RESOURCES/"

COPY src/WGD ${RESOURCES}WGD

COPY src/RdTest ${RESOURCES}RdTest

COPY src/svtk ${RESOURCES}svtk

COPY src/sv-pipeline ${RESOURCES}sv-pipeline

COPY src/svtest ${RESOURCES}svtest

COPY src/svqc ${RESOURCES}svqc
