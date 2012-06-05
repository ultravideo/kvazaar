/*! \file bitstream.h
    \brief Bitstream related functions
    \author Marko Viitanen
    \date 2012-05
    
    This file has all bitstream headers
*/
#ifndef _BITSTREAM_H
#define _BITSTREAM_H
 
 
typedef struct
{
    uint32_t data[32];
    uint8_t  cur_byte;
    uint8_t  cur_bit; 
    FILE*    output;
    uint8_t* buffer;
    uint32_t buffer_pos;
} bitstream;

typedef struct
{
  uint8_t len;
  uint32_t value;
}bitTable;

extern bitTable *exp_table;
 
void bitstream_alloc(bitstream* stream, uint32_t alloc); 
void bitstream_init(bitstream* stream); 
void bitstream_put(bitstream* stream, uint32_t data, uint8_t bits); 
/*
void bitstream_put_ue(bitstream* stream, uint32_t data);
void bitstream_put_se(bitstream* stream, uint32_t data);
*/
#define bitstream_put_ue(stream, data) { bitstream_put(stream,exp_table[data].value,exp_table[data].len); }
#define bitstream_put_se(stream, data) { uint32_t index=(data<=0)?2*(uint32_t)(-data):2*(uint32_t)(data)-1;    \
                                         bitstream_put(stream,exp_table[index].value,exp_table[index].len); }

void bitstream_align(bitstream* stream); 
void bitstream_flush(bitstream* stream);
void init_exp_golomb(uint32_t len);

#ifdef _DEBUG
static int WRITE_VALUE = 0;
#define WRITE_U(stream, data, bits, name) { printf("%8d %40s u(%d) : %d\n",WRITE_VALUE, name,bits,data); bitstream_put(stream,data,bits); WRITE_VALUE++;}
#define WRITE_UE(stream, data, name) { printf("%8d %40s ue(v) : %d\n",WRITE_VALUE, name,data); bitstream_put_ue(stream,data); WRITE_VALUE++;}
#define WRITE_SE(stream, data, name) { printf("%8d %40s se(v) : %d\n",WRITE_VALUE, name,data); bitstream_put_se(stream,data); WRITE_VALUE++;}
#else
#define WRITE_U(stream, data, bits, name) { bitstream_put(stream,data,bits); }
#define WRITE_UE(stream, data, name) { bitstream_put_ue(stream,data); }
#define WRITE_SE(stream, data, name) { bitstream_put_se(stream,data); }
#endif

 
#endif