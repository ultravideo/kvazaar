#ifndef KVAZAAR_VERSION_H
#define KVAZAAR_VERSION_H

#define KVA_VERSION_INT(a, b, c) (a<<16 | b<<8 | c)
#define KVA_VERSION_DOT(a, b, c) a ##.## b ##.## c
#define KVA_VERSION(a, b, c) KVA_VERSION_DOT(a, b, c)

#define KVAZAAR_VERSION_MAJOR  0
#define KVAZAAR_VERSION_MINOR  2
#define KVAZAAR_VERSION_MICRO  4

#define KVAZAAR_VERSION_INT   KVA_VERSION_INT(KVAZAAR_VERSION_MAJOR, \
                                               KVAZAAR_VERSION_MINOR, \
                                               KVAZAAR_VERSION_MICRO)
#define KVAZAAR_VERSION       KVA_VERSION(KVAZAAR_VERSION_MAJOR,     \
                                           KVAZAAR_VERSION_MINOR,     \
                                           KVAZAAR_VERSION_MICRO)
#define KVAZAAR_BUILD         KVAZAAR_VERSION_INT

#define KVA_STRING(s) # s

#define KVAZAAR_VERSION_STRING KVA_STRING(KVAZAAR_VERSION)

#endif /* KVAZAAR_VERSION_H */
