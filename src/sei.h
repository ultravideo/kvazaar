#ifndef UDU_SEI_H_
#define UDU_SEI_H_
/*****************************************************************************
 * This file is part of Kvazaar HEVC encoder.
 *
 * Copyright (c) 2022, Tampere University, ITU/ISO/IEC, project contributors
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

/**
 * \file
 * Supplemental Enhancement Information (SEI)
 */

#include "global.h" // IWYU pragma: keep

#define SEI_PAYLOAD_TYPE_PIC_TIMING 1
#define SEI_PAYLOAD_TYPE_USER_DATA_UNREGISTERED 5
#define SEI_PAYLOAD_TYPE_DECODED_PICTURE_HASH 132

// Flag value used for length / value extension.
#define FF_BYTE 0xFF

// UUID for the encoder info user defined unregistered message
// random uuid_iso_iec_11578 generated with www.famkruithof.net/uuid/uuidgen
static const uint8_t encoder_info_uuid[16] = {
    0x32, 0xfe, 0x46, 0x6c, 0x98, 0x41, 0x42, 0x69,
    0xae, 0x35, 0x6a, 0x91, 0x54, 0x9e, 0xf3, 0xf1};

#endif // UDU_SEI_H_
