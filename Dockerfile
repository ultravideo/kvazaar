FROM ubuntu:15.10
    ENV REQUIRED_PACKAGES automake autoconf libtool m4 build-essential git yasm
    RUN apt-get update
    RUN apt-get install -y $REQUIRED_PACKAGES ; \
    git clone --depth=1 git://github.com/ultravideo/kvazaar.git; \
        cd kvazaar; \
        ./autogen.sh; \
        ./configure --disable-shared;\
        make;\
    AUTOINSTALLED_PACKAGES=`apt-mark showauto`;\
    apt-get remove --purge -y $REQUIRED_PACKAGES $AUTOINSTALLED_PACKAGES; \
        apt-get autoremove -y; \
        apt-get clean all; \
        rm -rf /var/lib/apt/lists/*
    COPY src/kvazaar /
ENTRYPOINT ["/kvazaar"]
CMD ["--help"]
