FROM ubuntu:15.10
    RUN apt-get update
    RUN apt-get install -y \
        automake \
        autoconf \
        libtool \
        m4 \
        build-essential \
        git \
        yasm; \
    git clone --depth=1 git://github.com/ultravideo/kvazaar.git; \
        cd kvazaar; \
        ./autogen.sh; \
        ./configure --disable-shared;\
        make;\
    apt-get remove --purge -y automake autoconf libtool m4 build-essential git yasm `apt-mark showauto`; \
        apt-get autoremove -y; \
        apt-get clean all; \
        rm -rf /var/lib/apt/lists/*
    COPY src/kvazaar /
ENTRYPOINT ["/kvazaar"]
CMD ["--help"]
