# Continuous integration docker image
# Contains the bare files needed to build and test the application
FROM openjdk:17-bullseye as build
ENV BUILD_MODE 2
RUN apt update && apt install -y sudo

WORKDIR /opt/mog
COPY ./pom.xml ./pom.xml
COPY ./scripts ./scripts
RUN bash ./scripts/bootstrap.sh \
	&& mvn dependency:go-offline -B

COPY ./src ./src
RUN bash ./scripts/build.sh \
	&& cp ./target/metaomgraph4*.jar ./mog-applet.jar


# Live docker image
# 
#FROM openjdk:17-bullseye as live

# Graphics stuff
RUN apt update \
	&& apt install -y libxext6 libxrender1 libxtst6
#WORKDIR /opt/mog
#COPY --from=build /opt/mog/mog-applet.jar .
ENTRYPOINT ["java", "-jar", "./mog-applet.jar"]
