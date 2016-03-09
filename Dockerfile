FROM ubuntu:15.10
    RUN apt-get update
    RUN apt-get install -y \
        automake \
        autoconf \
        libtool \
        m4 \
        build-essential \
        git \
        yasm
    RUN git clone --depth=1 git://github.com/ultravideo/kvazaar.git
    RUN  cd kvazaar; \
        ./autogen.sh; \
        ./configure --disable-shared;\
        make
    COPY src/kvazaar /
    RUN apt-get purge -y automake autoconf libtool m4 build-essential git yasm
    RUN apt-get autoremove -y
    RUN apt-get clean all
    RUN rm -rf /var/lib/apt/lists/*
ENTRYPOINT ["/kvazaar"]
CMD ["--help"]
