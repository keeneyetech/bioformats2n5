# Development Dockerfile for bioformats2n5
# -----------------------------------------

# To install the built distribution into other runtimes
# pass a build argument, e.g.:
#
#   docker build --build-arg IMAGE=openjdk:9 ...
#

# Similarly, the BUILD_IMAGE argument can be overwritten
# but this is generally not needed.
ARG BUILD_IMAGE=gradle:6.2.1-jdk8

#
# Build phase: Use the gradle image for building.
#
FROM ${BUILD_IMAGE} as build
USER root
RUN apt-get update -qq && DEBIAN_FRONTEND=noninteractive apt-get install -y -qq zeroc-ice-all-runtime libblosc1
RUN mkdir /bioformats2n5 && chown 1000:1000 /bioformats2n5

# Build all
USER 1000

COPY --chown=1000:1000 . /bioformats2n5
WORKDIR /bioformats2n5
RUN gradle build
RUN cd build/distributions && rm bioformats2n5*tar && unzip bioformats2n5*zip && rm -rf bioformats2n5*zip
USER root
RUN mv /bioformats2n5/build/distributions/bioformats2n5* /opt/bioformats2n5
USER 1000
WORKDIR /opt/bioformats2n5
ENTRYPOINT ["/opt/bioformats2n5/bin/bioformats2n5"]
