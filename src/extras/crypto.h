
#include <stdio.h>
#include <math.h>
#define AESEncryptionStreamMode      1
#ifdef __cplusplus
extern "C" {
#endif
    typedef void* Crypto_Handle;
    Crypto_Handle InitC();
    void DecryptC(Crypto_Handle hdl, const unsigned char *in_stream, int size_bits, unsigned char  *out_stream);
#if AESEncryptionStreamMode
    unsigned int ff_get_key (Crypto_Handle *hdl, int nb_bits);
#endif
    void DeleteCryptoC(Crypto_Handle hdl);

#ifdef __cplusplus
}
#endif

