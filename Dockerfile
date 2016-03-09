# A simple Dockerfile for building Kvazaar from the git repository
# Example build command when in this directory: docker build -t kvazaar .

# Use Ubuntu 15.10 as a base for now, it's around 136MB
FROM ubuntu:15.10

    # List of needed packages to be able to build kvazaar with autotools
    ENV REQUIRED_PACKAGES automake autoconf libtool m4 build-essential git yasm
    
    # Running update causes extra data to be stored in the image
    # We don't need to update right now.
    #RUN apt-get update
    
    # Run all the commands in one RUN so we don't have any extra history
    # data in the image.
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
    # Because we build only the static binary, copy that to the root
    COPY src/kvazaar /
ENTRYPOINT ["/kvazaar"]
CMD ["--help"]
