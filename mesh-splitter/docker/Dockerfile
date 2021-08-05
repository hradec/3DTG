FROM ubuntu:20.04

RUN apt update
RUN apt install -y build-essential
RUN apt install -y sudo nano libmpfr-dev patchelf
RUN apt install -y assimp-utils libassimp5 libassimp-dev python3-pyassimp
RUN DEBIAN_FRONTEND="noninteractive" apt install -y tzdata
RUN apt install -y gdb
