/*****************************************************************************
 * This file is part of Kvazaar HEVC encoder.
 *
 * Copyright (c) 2023, Tampere University, ITU/ISO/IEC, project contributors
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above copyright notice, this
 *   list of conditions and the following disclaimer in the documentation and/or
 *   other materials provided with the distribution.
 *
 * * Neither the name of the Tampere University or ITU/ISO/IEC nor the names of its
 *   contributors may be used to endorse or promote products derived from
 *   this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * INCLUDING NEGLIGENCE OR OTHERWISE ARISING IN ANY WAY OUT OF THE USE OF THIS
 ****************************************************************************/

#include <extras/crypto.h>

#ifndef KVZ_SEL_ENCRYPTION
int kvz_make_vs_ignore_crypto_not_having_symbols = 0;
#else

#ifdef KVZ_ENABLE_BUILTIN_CRYPTO_SIMD
#define KVZ_CRYPTO_x86_SIMD 1
#endif

//#define KVZ_USE_CRYPTOPP 1

#if KVZ_CRYPTO_x86_SIMD == 1
#include <wmmintrin.h>
#endif

// Rcon table for AES.
const uint32_t aes_rcon[11] = {
    0x00000000, 0x01000000, 0x02000000,
    0x04000000, 0x08000000, 0x10000000,
    0x20000000, 0x40000000, 0x80000000,
    0x1b000000, 0x36000000
};

// AES S-Box data
const uint8_t aes_sbox[256] = {
    0x63, 0x7C, 0x77, 0x7B, 0xF2, 0x6B, 0x6F, 0xC5, 0x30, 0x01, 0x67, 0x2B, 0xFE, 0xD7, 0xAB, 0x76,
    0xCA, 0x82, 0xC9, 0x7D, 0xFA, 0x59, 0x47, 0xF0, 0xAD, 0xD4, 0xA2, 0xAF, 0x9C, 0xA4, 0x72, 0xC0,
    0xB7, 0xFD, 0x93, 0x26, 0x36, 0x3F, 0xF7, 0xCC, 0x34, 0xA5, 0xE5, 0xF1, 0x71, 0xD8, 0x31, 0x15,
    0x04, 0xC7, 0x23, 0xC3, 0x18, 0x96, 0x05, 0x9A, 0x07, 0x12, 0x80, 0xE2, 0xEB, 0x27, 0xB2, 0x75,
    0x09, 0x83, 0x2C, 0x1A, 0x1B, 0x6E, 0x5A, 0xA0, 0x52, 0x3B, 0xD6, 0xB3, 0x29, 0xE3, 0x2F, 0x84,
    0x53, 0xD1, 0x00, 0xED, 0x20, 0xFC, 0xB1, 0x5B, 0x6A, 0xCB, 0xBE, 0x39, 0x4A, 0x4C, 0x58, 0xCF,
    0xD0, 0xEF, 0xAA, 0xFB, 0x43, 0x4D, 0x33, 0x85, 0x45, 0xF9, 0x02, 0x7F, 0x50, 0x3C, 0x9F, 0xA8,
    0x51, 0xA3, 0x40, 0x8F, 0x92, 0x9D, 0x38, 0xF5, 0xBC, 0xB6, 0xDA, 0x21, 0x10, 0xFF, 0xF3, 0xD2,
    0xCD, 0x0C, 0x13, 0xEC, 0x5F, 0x97, 0x44, 0x17, 0xC4, 0xA7, 0x7E, 0x3D, 0x64, 0x5D, 0x19, 0x73,
    0x60, 0x81, 0x4F, 0xDC, 0x22, 0x2A, 0x90, 0x88, 0x46, 0xEE, 0xB8, 0x14, 0xDE, 0x5E, 0x0B, 0xDB,
    0xE0, 0x32, 0x3A, 0x0A, 0x49, 0x06, 0x24, 0x5C, 0xC2, 0xD3, 0xAC, 0x62, 0x91, 0x95, 0xE4, 0x79,
    0xE7, 0xC8, 0x37, 0x6D, 0x8D, 0xD5, 0x4E, 0xA9, 0x6C, 0x56, 0xF4, 0xEA, 0x65, 0x7A, 0xAE, 0x08,
    0xBA, 0x78, 0x25, 0x2E, 0x1C, 0xA6, 0xB4, 0xC6, 0xE8, 0xDD, 0x74, 0x1F, 0x4B, 0xBD, 0x8B, 0x8A,
    0x70, 0x3E, 0xB5, 0x66, 0x48, 0x03, 0xF6, 0x0E, 0x61, 0x35, 0x57, 0xB9, 0x86, 0xC1, 0x1D, 0x9E,
    0xE1, 0xF8, 0x98, 0x11, 0x69, 0xD9, 0x8E, 0x94, 0x9B, 0x1E, 0x87, 0xE9, 0xCE, 0x55, 0x28, 0xDF,
    0x8C, 0xA1, 0x89, 0x0D, 0xBF, 0xE6, 0x42, 0x68, 0x41, 0x99, 0x2D, 0x0F, 0xB0, 0x54, 0xBB, 0x16
};


static void aes_sub_bytes(uint8_t* state) {
  for (int i = 0; i < 16; i++) {
    state[i] = aes_sbox[state[i]];
  }
}


static void aes_shift_rows(uint8_t* state) {
  uint8_t temp[16];
  memcpy(temp, state, 16);

  state[1] = temp[5]; state[5] = temp[9];  state[9] = temp[13]; state[13] = temp[1];
  state[2] = temp[10]; state[6] = temp[14]; state[10] = temp[2]; state[14] = temp[6];
  state[3] = temp[15]; state[7] = temp[3]; state[11] = temp[7]; state[15] = temp[11];
}


/*
This function performs multiplication in the Galois field GF(2^8)
| 2 3 1 1 |
| 1 2 3 1 |
| 1 1 2 3 |
| 3 1 1 2 |
*/
static uint8_t aes_gmul(uint8_t a, uint8_t b) {
  uint8_t p = 0;
  const uint8_t high_bit_mask = 0x80;
  uint8_t high_bit = 0;
  const uint8_t modulo = 0x1B; /* x^8 + x^4 + x^3 + x + 1 */

  for (int i = 0; i < 8; i++) {
    if (b & 1) p ^= a;
    high_bit = a & high_bit_mask;
    a <<= 1;
    if (high_bit) a ^= modulo;
    b >>= 1;
  }
  return p;
}

static void aes_mix_columns(uint8_t* state) {
  uint8_t tmp[16];
  memcpy(tmp, state, 16);
  for (int i = 0; i < 4; ++i) {
    state[i * 4 + 0] = aes_gmul(tmp[i * 4 + 0], 0x02) ^ aes_gmul(tmp[i * 4 + 1], 0x03) ^ tmp[i * 4 + 2] ^ tmp[i * 4 + 3];
    state[i * 4 + 1] = tmp[i * 4 + 0] ^ aes_gmul(tmp[i * 4 + 1], 0x02) ^ aes_gmul(tmp[i * 4 + 2], 0x03) ^ tmp[i * 4 + 3];
    state[i * 4 + 2] = tmp[i * 4 + 0] ^ tmp[i * 4 + 1] ^ aes_gmul(tmp[i * 4 + 2], 0x02) ^ aes_gmul(tmp[i * 4 + 3], 0x03);
    state[i * 4 + 3] = aes_gmul(tmp[i * 4 + 0], 0x03) ^ tmp[i * 4 + 1] ^ tmp[i * 4 + 2] ^ aes_gmul(tmp[i * 4 + 3], 0x02);
  }
}


static void aes_add_roundkey(uint8_t* state, uint8_t* roundKey) {
  for (int i = 0; i < 16; ++i) {
    state[i] ^= roundKey[i];
  }
}


// SubWord for uint32_t, replaces each byte with the corresponding byte from the sBox
static uint32_t aes_sub_word(uint32_t word) {
  uint8_t byte[4] = { (uint8_t)(word >> 24), (uint8_t)(word >> 16), (uint8_t)(word >> 8), (uint8_t)(word) };
  word = (aes_sbox[byte[0]] << 24) | (aes_sbox[byte[1]] << 16) | (aes_sbox[byte[2]] << 8) | (aes_sbox[byte[3]]);
  return word;
}
static uint32_t aes_rot_word(uint32_t word) {
#ifndef __BIG_ENDIAN__
  return (word << 24) | (word >> 8);
#else
  return (word << 8) | (word >> 24);
#endif
}

// The key expansion function, generates 11 round keys from the cipher key
static uint32_t* aes_key_expansion(const uint8_t* cipherKey, uint32_t* w) {

  for (int i = 0; i < 4; ++i) {
#ifndef __BIG_ENDIAN__
    w[i] = (cipherKey[i * 4 + 3] << 24) | (cipherKey[i * 4 + 2] << 16) | (cipherKey[i * 4 + 1] << 8) | (cipherKey[i * 4 + 0]);
#else
    w[i] = (cipherKey[i * 4] << 24) | (cipherKey[i * 4 + 1] << 16) | (cipherKey[i * 4 + 2] << 8) | (cipherKey[i * 4 + 3]);
#endif
  }

  for (int i = 4; i < 44; ++i) {
    uint32_t temp = w[i - 1];
    if (i % 4 == 0) {
      temp = aes_sub_word(aes_rot_word(temp));
#ifndef __BIG_ENDIAN__
      temp = temp ^ (aes_rcon[i / 4] >> 24);
#else
      temp = temp ^ (aes_rcon[i / 4]);
#endif
    }
    w[i] = w[i - 4] ^ temp;
  }
  return w;
}

static void aes128_encrypt_with_roundkeys(const uint8_t* plaintext, uint8_t* ciphertext, uint32_t* roundKeys) {
  uint8_t state[16];
  memcpy(state, plaintext, 16);

  // Initial AddRoundKey step
  aes_add_roundkey(state, (uint8_t*)&roundKeys[0]);

  // 9 rounds of SubBytes, ShiftRows, MixColumns, and AddRoundKey
  for (int round = 1; round <= 9; ++round) {
#if KVZ_CRYPTO_x86_SIMD == 1
    * (__m128i*)state = _mm_aesenc_si128(_mm_loadu_si128((__m128i*)state), _mm_loadu_si128((__m128i*)&roundKeys[round * 4]));
#else
    aes_sub_bytes(state);
    aes_shift_rows(state);
    aes_mix_columns(state);
    aes_add_roundkey(state, (uint8_t*)&roundKeys[round * 4]);
#endif
  }
#if KVZ_CRYPTO_x86_SIMD == 1
  *(__m128i*)state = _mm_aesenclast_si128(_mm_loadu_si128((__m128i*)state), _mm_loadu_si128((__m128i*) & roundKeys[10 * 4]));
#else
  // Final round (no MixColumns)
  aes_sub_bytes(state);
  aes_shift_rows(state);
  aes_add_roundkey(state, (uint8_t*)&roundKeys[10 * 4]);
#endif
  memcpy(ciphertext, state, 16);
}

static void aes128_encrypt_CFB_with_roundkeys(const uint8_t* plaintext, uint8_t* ciphertext, uint8_t* IV, uint32_t* roundKeys) {
  // AES-128 encryption of the IV
  uint8_t encrypted_IV[16];
  aes128_encrypt_with_roundkeys(IV, encrypted_IV,roundKeys);

  // XOR plaintext with encrypted IV
  for (int i = 0; i < 16; ++i) {
    IV[i] = plaintext[i] ^ encrypted_IV[i];
  }
  memcpy(ciphertext, IV, 16);
}





#ifdef KVZ_USE_CRYPTOPP

#include <cryptopp/aes.h>
#include <cryptopp/modes.h>
#include <cryptopp/osrng.h>

#if AESEncryptionStreamMode
  typedef CryptoPP::CFB_Mode<CryptoPP::AES>::Encryption cipher_t;
#else
  typedef CryptoPP::CFB_Mode<CryptoPP::AES>::Decryption cipher_t;
#endif
#endif
struct crypto_handle_t {
#ifdef KVZ_USE_CRYPTOPP
  cipher_t *cipher;
  unsigned char key[CryptoPP::AES::DEFAULT_KEYLENGTH];
  unsigned char iv[CryptoPP::AES::BLOCKSIZE];
  unsigned char out_stream_counter[CryptoPP::AES::BLOCKSIZE];
  unsigned char counter[CryptoPP::AES::BLOCKSIZE];
#else
  uint8_t key[16];
  uint8_t iv[16];
  uint8_t out_stream_counter[16];
  uint8_t counter[16];
  uint8_t kvz_iv[16];
  uint32_t kvz_roundkeys[44];
#endif
  int couter_avail;
  int counter_index;
  int counter_index_pos;
};


static uint8_t default_IV[16] = {201, 75, 219, 152, 6, 245, 237, 107, 179, 194, 81, 29, 66, 98, 198, 0};
static uint8_t default_key[16] = {16, 213, 27, 56, 255, 127, 242, 112, 97, 126, 197, 204, 25, 59, 38, 30};


crypto_handle_t* kvz_crypto_create(const kvz_config *cfg)
{
  crypto_handle_t* hdl = (crypto_handle_t*)calloc(1, sizeof(crypto_handle_t));

  uint8_t *key;
  if(cfg->optional_key!=NULL)
    key = cfg->optional_key;
  else
    key = default_key;

  for (int i = 0; i < 16; i++) {
    hdl->iv [i]     = default_IV[i];
    hdl->counter[i] = (i<11)? default_IV[5+i] : key[i-11];
    hdl->key[i]     = key[i];
  }
    
#ifdef KVZ_USE_CRYPTOPP
  hdl->cipher = new cipher_t(hdl->key, CryptoPP::AES::DEFAULT_KEYLENGTH, hdl->iv);
#else
  memcpy(hdl->kvz_iv, hdl->iv, 16); 
  aes_key_expansion(hdl->key, hdl->kvz_roundkeys);
#endif
  hdl->couter_avail      = 0;
  hdl->counter_index     = 0;
  hdl->counter_index_pos = 0;

  return hdl;
}

void kvz_crypto_delete(crypto_handle_t **hdl)
{
#ifdef KVZ_USE_CRYPTOPP
  if (*hdl) {
    delete (*hdl)->cipher;
    (*hdl)->cipher = NULL;
  }
#endif
  FREE_POINTER(*hdl);
}

void kvz_crypto_decrypt(crypto_handle_t* hdl,
                        const uint8_t *in_stream,
                        int size_bits,
                        uint8_t *out_stream)
{
  int num_bytes = ceil((double)size_bits/8);
#ifdef KVZ_USE_CRYPTOPP
  hdl->cipher->ProcessData(out_stream, in_stream, num_bytes);
  if (size_bits & 7) {
    hdl->cipher->SetKeyWithIV(hdl->key, CryptoPP::AES::DEFAULT_KEYLENGTH, hdl->iv);
  }
#else
  // Hopefully handles the same case as cryptopp
  uint32_t bytes_left = num_bytes;
  uint32_t buffer_pos = 0;
  while (bytes_left >= 16) {
    aes128_encrypt_CFB_with_roundkeys(in_stream, out_stream, hdl->kvz_iv, hdl->kvz_roundkeys);
    buffer_pos += 16;
  }
  if (bytes_left > 0) {
    uint8_t block_in[16] = { 0 };
    uint8_t block_out[16] = { 0 };
    memcpy(block_in, &in_stream[buffer_pos], bytes_left);
    aes128_encrypt_CFB_with_roundkeys(block_in, block_out, hdl->kvz_iv, hdl->kvz_roundkeys);
    memcpy(&out_stream[buffer_pos], block_out, bytes_left);
  }
  if (size_bits & 7) {
    memcpy(hdl->kvz_iv, hdl->iv, 16);
  }
#endif
}
#if AESEncryptionStreamMode
static void increment_counter(unsigned char *counter)
{
  counter[0]++;
}

static void decrypt_counter(crypto_handle_t *hdl)
{
#ifdef KVZ_USE_CRYPTOPP
  hdl->cipher->ProcessData(hdl->out_stream_counter, hdl->counter, 16);
#else
  aes128_encrypt_CFB_with_roundkeys(hdl->counter, hdl->out_stream_counter,hdl->kvz_iv, hdl->kvz_roundkeys);
#endif
  hdl->couter_avail      = 128;
  hdl->counter_index     = 15;
  hdl->counter_index_pos = 8;
  increment_counter(hdl->counter);
}

unsigned kvz_crypto_get_key(crypto_handle_t *hdl, int nb_bits)
{
  unsigned key = 0;
  if (nb_bits > 32) {
      fprintf(stderr, "The generator cannot generate %d bits (max 32 bits)\n", nb_bits);
      return 0;
  }
  if (nb_bits == 0) return 0;

  if (!hdl->couter_avail) {
    decrypt_counter(hdl);
  }

  if(hdl->couter_avail >= nb_bits) {
      hdl->couter_avail -= nb_bits;
  } else {
      hdl->couter_avail = 0;
  }

  int nb = 0;
  while (nb_bits) {
    if (nb_bits >= hdl->counter_index_pos) {
      nb = hdl->counter_index_pos;
    } else {
      nb = nb_bits;
    }

    key <<= nb;
    key += hdl->out_stream_counter[hdl->counter_index] & ((1 << nb) - 1);
    hdl->out_stream_counter[hdl->counter_index] >>= nb;
    nb_bits -= nb;

    if (hdl->counter_index && nb == hdl->counter_index_pos) {
      hdl->counter_index--;
      hdl->counter_index_pos = 8;
    } else {
      hdl->counter_index_pos -= nb;
      if (nb_bits) {
        decrypt_counter(hdl);
        hdl->couter_avail -=  nb_bits;
      }
    }
  }
  return key;
}
#endif // AESEncryptionStreamMode

#endif // KVZ_SEL_ENCRYPTION
