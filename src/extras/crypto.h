#ifndef CRYPTO_H_
#define CRYPTO_H_

#include "global.h"

#ifdef KVZ_SEL_ENCRYPTION
#define STUBBED extern
#else
#define STUBBED static
#endif

#include <stdio.h>
#include <math.h>
#define AESEncryptionStreamMode      1
#ifdef __cplusplus
extern "C" {
#endif
    typedef void* Crypto_Handle;

    STUBBED Crypto_Handle InitC();
    STUBBED void DecryptC(Crypto_Handle hdl, const unsigned char *in_stream, int size_bits, unsigned char  *out_stream);
#if AESEncryptionStreamMode
    STUBBED unsigned int ff_get_key(Crypto_Handle *hdl, int nb_bits);
#endif
    STUBBED void DeleteCryptoC(Crypto_Handle hdl);

#ifdef __cplusplus
}
#endif


#ifndef KVZ_SEL_ENCRYPTION
// Provide static stubs to allow linking without libcryptopp and allows us to
// avoid sprinkling ifdefs everywhere and having a bunch of code that's not
// compiled during normal development.
// Provide them in the header so we can avoid compiling the cpp file, which
// means we don't need a C++ compiler when crypto is not enabled.

#include <assert.h>

static INLINE Crypto_Handle InitC()
{
  // Stub.
  assert(0);
  return 0;
}

static INLINE void DecryptC(Crypto_Handle hdl, const unsigned char *in_stream,
              int size_bits, unsigned char  *out_stream)
{
  // Stub.
  assert(0);
}

#if AESEncryptionStreamMode
static INLINE unsigned int ff_get_key(Crypto_Handle *hdl, int nb_bits)
{
  // Stub.
  assert(0);
  return 0;
}
#endif

static INLINE void DeleteCryptoC(Crypto_Handle hdl)
{
  // Stub.
  assert(0);
}

#endif // KVZ_SEL_ENCRYPTION

#endif // CRYPTO_H_
