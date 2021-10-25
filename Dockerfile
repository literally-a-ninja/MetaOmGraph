# Continuous integration docker image
#
# Contains the bare files needed to build and test the application
FROM openjdk:17-bullseye as ci
RUN apt update
RUN apt install -y maven wget unzip
WORKDIR /usr/share/mog

COPY pom.xml pom.xml
COPY src/ src/
COPY scripts/ scripts/

RUN bash scripts/build.sh
RUN bash scripts/test.sh

ENTRYPOINT java -jar target/metaomgraph4-1.8.2beta-jar-with-dependencies.jar
