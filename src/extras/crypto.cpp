#include <extras/crypto.h>
#include <cryptopp/aes.h>
#include <cryptopp/modes.h>
#include <cryptopp/osrng.h>
typedef struct AESDecoder {
#if AESEncryptionStreamMode
        CryptoPP::CFB_Mode<CryptoPP::AES>::Encryption *CFBdec;
#else
    CryptoPP::CFB_Mode<CryptoPP::AES>::Decryption *CFBdec;
#endif

    byte key[CryptoPP::AES::DEFAULT_KEYLENGTH], iv[CryptoPP::AES::BLOCKSIZE], out_stream_counter[CryptoPP::AES::BLOCKSIZE], counter[CryptoPP::AES::BLOCKSIZE];
    int couter_avail, counter_index, counter_index_pos;
} AESDecoder;


AESDecoder* Init() {
    int init_val[32] = {201, 75, 219, 152, 6, 245, 237, 107, 179, 194, 81, 29, 66, 98, 198, 0, 16, 213, 27, 56, 255, 127, 242, 112, 97, 126, 197, 204, 25, 59, 38, 30};
    AESDecoder * AESdecoder = (AESDecoder *)malloc(sizeof(AESDecoder));
    for(int i=0;i<16; i++) {
        AESdecoder->iv [i]     = init_val[i];
        AESdecoder->counter[i] = init_val[5+i];
        AESdecoder->key[i]     = init_val[i+16];
    }
#if AESEncryptionStreamMode
    AESdecoder->CFBdec = new CryptoPP::CFB_Mode<CryptoPP::AES >::Encryption(AESdecoder->key, CryptoPP::AES::DEFAULT_KEYLENGTH, AESdecoder->iv);
#else
    AESdecoder->CFBdec = new CryptoPP::CFB_Mode<CryptoPP::AES >::Decryption(AESdecoder->key, CryptoPP::AES::DEFAULT_KEYLENGTH, AESdecoder->iv);
#endif
    AESdecoder->couter_avail      = 0;
    AESdecoder->counter_index     = 0;
    AESdecoder->counter_index_pos = 0;
    return AESdecoder;
}

void DeleteCrypto(AESDecoder * AESdecoder) {
    if(AESdecoder)
        free(AESdecoder);
}

void Decrypt(AESDecoder *AESdecoder, const unsigned char *in_stream, int size_bits, unsigned char  *out_stream) {
    int nb_bytes = ceil((double)size_bits/8);
    AESdecoder->CFBdec->ProcessData(out_stream, in_stream, nb_bytes);
    if(size_bits&7)
        AESdecoder->CFBdec->SetKeyWithIV(AESdecoder->key, CryptoPP::AES::DEFAULT_KEYLENGTH, AESdecoder->iv);
    
}
void Incr_counter (unsigned char *counter) {
    counter[0]++;
}

#if AESEncryptionStreamMode
void Decrypt_counter(AESDecoder * AESdecoder) {
    AESdecoder->CFBdec->ProcessData(AESdecoder->out_stream_counter, AESdecoder->counter, 16);
    AESdecoder->couter_avail      = 128;
    AESdecoder->counter_index     = 15;
    AESdecoder->counter_index_pos = 8;
    Incr_counter(AESdecoder->counter);
}
#endif

#if AESEncryptionStreamMode
unsigned int get_key (AESDecoder * AESdecoder, int nb_bits) {
    unsigned int key_ = 0;
    if(nb_bits > 32) {
        printf("The Generator can not generate more than 32 bit %d \n", nb_bits);
        return 0;
    }
    if( !nb_bits )
        return 0;
    if(!AESdecoder->couter_avail)
        Decrypt_counter(AESdecoder);

    if(AESdecoder->couter_avail >= nb_bits)
        AESdecoder->couter_avail -= nb_bits;
    else
        AESdecoder->couter_avail = 0;
    int nb = 0;
    while( nb_bits ) {
        if( nb_bits >= AESdecoder->counter_index_pos )
            nb = AESdecoder->counter_index_pos;
        else
            nb = nb_bits;
        key_ <<= nb;
        key_ += (AESdecoder->out_stream_counter[AESdecoder->counter_index] & ((1<<nb)-1));
        AESdecoder->out_stream_counter[AESdecoder->counter_index] >>= nb;
        nb_bits -= nb;

        if(AESdecoder->counter_index && nb == AESdecoder->counter_index_pos ) {
            AESdecoder->counter_index--;
            AESdecoder->counter_index_pos = 8;
        } else {
            AESdecoder->counter_index_pos -= nb;
            if(nb_bits) {
                Decrypt_counter(AESdecoder);
                AESdecoder->couter_avail -=  nb_bits;
            }
        }
    }
    return key_;
}
#endif

Crypto_Handle InitC(){
    AESDecoder* AESdecoder = Init();
    return AESdecoder;
}
#if AESEncryptionStreamMode
unsigned int ff_get_key (Crypto_Handle *hdl, int nb_bits) {
    return get_key ((AESDecoder*)*hdl, nb_bits);
}
#endif
void DecryptC(Crypto_Handle hdl, const unsigned char *in_stream, int size_bits, unsigned char  *out_stream) {
    Decrypt((AESDecoder*)hdl, in_stream, size_bits, out_stream);
}

void DeleteCryptoC(Crypto_Handle hdl) {
    DeleteCryptoC(hdl);
}
