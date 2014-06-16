;/*****************************************************************************
;* This file is part of Kvazaar HEVC encoder.
;* 
;* Copyright (C) 2013-2014 Tampere University of Technology and others (see 
;* COPYING file).
;*
;* Kvazaar is free software: you can redistribute it and/or modify
;* it under the terms of the GNU General Public License version 2 as published
;* by the Free Software Foundation.
;*
;* Kvazaar is distributed in the hope that it will be useful,
;* but WITHOUT ANY WARRANTY; without even the implied warranty of
;* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;* GNU General Public License for more details.
;*
;* You should have received a copy of the GNU General Public License
;* along with Kvazaar.  If not, see <http://www.gnu.org/licenses/>.
;****************************************************************************/

%include "x86inc.asm"



SECTION .text

;KVZ_SAD_4X4
;Calculates SAD of the 16 consequtive bytes in memory
;r0 address of the first value(current)
;r1 address of the first value(reference)

cglobal sad_4x4, 2, 2, 2

vmovdqu m0, [r0]
vmovdqu m1, [r1]

vpsadbw m0, m1

vmovhlps m1, m0
vpaddw m0, m1

vmovd eax, m0

RET


;KVZ_SAD_4X4_STRIDE
;Calculates SAD of a 4x4 block inside a frame with stride
;r0 address of the first value(current)
;r1 address of the first value(reference)
;r2 stride

cglobal sad_4x4_stride, 3, 3, 2

vpinsrd m0, [r0], 0
add r0, r2
vpinsrd m0, [r0], 1
vpinsrd m0, [r0+r2], 2
vpinsrd m0, [r0+r2*2], 3

vpinsrd m1, [r1], 0
add r1, r2
vpinsrd m1, [r1], 1
vpinsrd m1, [r1+r2], 2
vpinsrd m1, [r1+r2*2], 3

vpsadbw m0, m1

vmovhlps m1, m0
vpaddw m0, m1

vmovd eax, m0

RET


;KVZ_SAD_8X8
;Calculates SAD of the 64 consequtive bytes in memory
;r0 address of the first value(current)
;r1 address of the first value(reference)

cglobal sad_8x8, 2, 2, 5

vpxor m0, m0

%rep 2

vmovdqu m1, [r0]
vmovdqu m3, [r0+16]
add r0, 32

vmovdqu m2, [r1]
vmovdqu m4, [r1+16]
add r1, 32

vpsadbw m1, m2
vpsadbw m3, m4

vpaddw m0, m1
vpaddw m0, m3

%endrep

vmovhlps m1, m0
vpaddw m0, m1

vmovd eax, m0

RET


;KVZ_SAD_8X8_STRIDE
;Calculates SAD of a 8x8 block inside a frame with stride
;r0 address of the first value(current)
;r1 address of the first value(reference)
;r2 stride

cglobal sad_8x8_stride, 3, 3, 5

vpxor m0, m0

vmovhpd m1, [r0]
add r0, r2
vmovlpd m1, [r0] 

vmovhpd m3, [r0+r2]
vmovlpd m3, [r0+r2*2] 
lea r0, [r0+r2*2]
add r0, r2

vmovhpd m2, [r1]
add r1, r2
vmovlpd m2, [r1] 

vmovhpd m4, [r1+r2]
vmovlpd m4, [r1+r2*2] 
lea r1, [r1+r2*2]
add r1, r2

vpsadbw m1, m2
vpsadbw m3, m4

vpaddw m0, m1
vpaddw m0, m3

vmovhpd m1, [r0]
add r0, r2
vmovlpd m1, [r0] 

vmovhpd m3, [r0+r2]
vmovlpd m3, [r0+r2*2] 
lea r0, [r0+r2*2]
add r0, r2

vmovhpd m2, [r1]
add r1, r2
vmovlpd m2, [r1] 

vmovhpd m4, [r1+r2]
vmovlpd m4, [r1+r2*2] 
lea r1, [r1+r2*2]
add r1, r2

vpsadbw m1, m2
vpsadbw m3, m4

vpaddw m0, m1
vpaddw m0, m3

vmovhlps m1, m0
vpaddw m0, m1

vmovd eax, m0

RET


;KVZ_SAD_16X16
;Calculates SAD of the 256 consequtive bytes in memory
;r0 address of the first value(current)
;r1 address of the first value(reference)

cglobal sad_16x16, 2, 2, 5

vpxor m4, m4

%rep 8

; Load 2 rows from rec_buf to m0 and m2
vmovdqu m0, [r0]
vmovdqu m2, [r0 + 16]
add r0, 32

; Load 2 rows from ref_buf to m1 and m3
vmovdqu m1, [r1]
vmovdqu m3, [r1 + 16]
add r1, 32

vpsadbw m0, m1
vpsadbw m2, m3

vpaddw m4, m0
vpaddw m4, m2

%endrep

vmovhlps m0, m4
vpaddw m4, m0

vmovd eax, m4

RET


;KVZ_SAD_16X16_STRIDE
;Calculates SAD of a 16x16 block inside a frame with stride
;r0 address of the first value(current)
;r1 address of the first value(reference)
;r2 stride

cglobal sad_16x16_stride, 3, 3, 5

vpxor m4, m4

%rep 8

; Load 2 rows from rec_buf to m0 and m2
vmovdqu m0, [r0]
vmovdqu m2, [r0 + r2]
lea r0, [r0 + r2*2]

; Load 2 rows from ref_buf to m1 and m3
vmovdqu m1, [r1]
vmovdqu m3, [r1 + r2]
lea r1, [r1 + r2*2]
 
vpsadbw m0, m1
vpsadbw m2, m3

vpaddw m4, m0
vpaddw m4, m2

%endrep

vmovhlps m0, m4
vpaddw m4, m0

vmovd eax, m4

RET


;KVZ_SATD_4X4
;Calculates SATD of the 16 consequtive bytes in memory
;r0 address of the first value(current)
;r1 address of the first value(reference)

cglobal satd_4x4, 2, 2, 6

vpmovzxbw m0, [r0]
vpmovzxbw m2, [r1]
vpsubw m0, m2

vpmovzxbw m1, [r0+8]
vpmovzxbw m3, [r1+8]
vpsubw m1, m3

;Horizontal phase
;rows 1-2
vphaddw m4, m0, m1
vphsubw m5, m0, m1

vphaddw m0, m4, m5
vphsubw m1, m4, m5

;Vertical phase
vphaddw m4, m0, m1
vphsubw m5, m0, m1

vphaddw m0, m4, m5
vphsubw m1, m4, m5

vpabsw m0, m0
vpabsw m1, m1

vpaddw m0, m1

vphaddw m0, m0
vphaddw m0, m0
vphaddw m0, m0

vpextrw eax, m0, 0

;Uncomment if transformed values not divided elsewhere
;add eax, 1
;shr eax, 1

RET


;KVZ_SATD_8X8_STRIDE
;Calculates SATD of a 8x8 block inside a frame with stride
;r0 address of the first value(reference)
;r1 address of the first value(current)
;r2 stride

%if ARCH_X86_64
    cglobal satd_8x8_stride, 4, 4, 16
%else
    cglobal satd_8x8_stride, 4, 4, 8
%endif

vpmovzxbw m0, [r0]
vpmovzxbw m7, [r2]
vpsubw m0, m7

vpmovzxbw m1, [r0+r1]
lea r0, [r0+r1*2]
vpmovzxbw m7, [r2+r3]
lea r2, [r2+r3*2]
vpsubw m1, m7

vpmovzxbw m2, [r0]
vpmovzxbw m7, [r2]
vpsubw m2, m7

vpmovzxbw m3, [r0+r1]
lea r0, [r0+r1*2]
vpmovzxbw m7, [r2+r3]
lea r2, [r2+r3*2]
vpsubw m3, m7

vpmovzxbw m4, [r0]
vpmovzxbw m7, [r2]
vpsubw m4, m7

vpmovzxbw m5, [r0+r1]
lea r0, [r0+r1*2]
vpmovzxbw m7, [r2+r3]
lea r2, [r2+r3*2]
vpsubw m5, m7

vpmovzxbw m6, [r0]
vpmovzxbw m7, [r2]
vpsubw m6, m7


%if ARCH_X86_64
    vpmovzxbw m7, [r0+r1]
    vpmovzxbw m8, [r2+r3]
    vpsubw m7, m8
%elif
    vpmovzxbw m7, [r2+r3]
    movdqu [esp-16], m7
    vpmovzxbw m7, [r0+r1]
    vpsubw m7, [esp-16]

    movdqu [esp-16], m4
    movdqu [esp-16*2], m5
    movdqu [esp-16*3], m6
    movdqu [esp-16*4], m7
    lea esp, [esp-16*4]
%endif

;Horizontal phaze

%if ARCH_X86_64
vphaddw m8, m0, m1
vphsubw m9, m0, m1

vphaddw m10, m2, m3
vphsubw m11, m2, m3

vphaddw m12, m4, m5
vphsubw m13, m4, m5

vphaddw m14, m6, m7 
vphsubw m15, m6, m7


vphaddw m0, m8, m9
vphsubw m1, m8, m9

vphaddw m2, m10, m11
vphsubw m3, m10, m11

vphaddw m4, m12, m13
vphsubw m5, m12, m13

vphaddw m6, m14, m15
vphsubw m7, m14, m15


vphaddw m8, m0, m1
vphsubw m9, m0, m1

vphaddw m10, m2, m3
vphsubw m11, m2, m3

vphaddw m12, m4, m5
vphsubw m13, m4, m5

vphaddw m14, m6, m7
vphsubw m15, m6, m7

%elif

vphaddw m4, m0, m1
vphsubw m5, m0, m1

vphaddw m6, m2, m3
vphsubw m7, m2, m3


vphaddw m0, m4, m5
vphsubw m1, m4, m5

vphaddw m2, m6, m7
vphsubw m3, m6, m7


vphaddw m4, m0, m1
vphsubw m5, m0, m1

vphaddw m6, m2, m3
vphsubw m7, m2, m3


movdqu m3, [esp]
movdqu m2, [esp+16]
movdqu m1, [esp+16*2]
movdqu m0, [esp+16*3]

movdqu [esp], m7
movdqu [esp+16*1], m6
movdqu [esp+16*2], m5
movdqu [esp+16*3], m4

vphaddw m4, m0, m1
vphsubw m5, m0, m1

vphaddw m6, m2, m3
vphsubw m7, m2, m3


vphaddw m0, m4, m5
vphsubw m1, m4, m5

vphaddw m2, m6, m7
vphsubw m3, m6, m7


vphaddw m4, m0, m1
vphsubw m5, m0, m1

vphaddw m6, m2, m3
vphsubw m7, m2, m3

%endif


;Vertical phase

%if ARCH_X86_64

vphaddw m0, m8, m9
vphsubw m1, m8, m9

vphaddw m2, m10, m11
vphsubw m3, m10, m11

vphaddw m4, m12, m13
vphsubw m5, m12, m13

vphaddw m6, m14, m15
vphsubw m7, m14, m15

vpmovzxwd m0, m0
vpmovzxwd m1, m1
vpmovzxwd m2, m2
vpmovzxwd m3, m3
vpmovzxwd m4, m4
vpmovzxwd m5, m5
vpmovzxwd m6, m6
vpmovzxwd m7, m7

vpaddd m8, m0, m2
vpaddd m9, m1, m3
vpsubd m10, m0, m2
vpsubd m11, m1, m3

vpaddd m12, m4, m6
vpaddd m13, m5, m7
vpsubd m14, m4, m6
vpsubd m15, m5, m7

vpaddd m0, m8, m12
vpaddd m1, m9, m13
vpaddd m2, m10, m14
vpaddd m3, m11, m15

vpsubd m4, m8, m12
vpsubd m5, m9, m13
vpsubd m6, m10, m14
vpsubd m7, m11, m15

%elif

vphaddw m0, m4, m5
vphsubw m1, m4, m5

vphaddw m2, m6, m7
vphsubw m3, m6, m7

vpmovzxwd m0, m0
vpmovzxwd m1, m1
vpmovzxwd m2, m2
vpmovzxwd m3, m3

vpaddd m4, m0, m2
vpaddd m5, m1, m3
vpsubd m6, m0, m2
vpsubd m7, m1, m3

movdqu m3, [esp]
movdqu m2, [esp+16]
movdqu m1, [esp+16*2]
movdqu m0, [esp+16*3]

movdqu [esp], m7
movdqu [esp+16*1], m6
movdqu [esp+16*2], m5
movdqu [esp+16*3], m4

vphaddw m0, m4, m5
vphsubw m1, m4, m5

vphaddw m2, m6, m7
vphsubw m3, m6, m7

vpmovzxwd m0, m0
vpmovzxwd m1, m1
vpmovzxwd m2, m2
vpmovzxwd m3, m3

vpaddd m4, m0, m2
vpaddd m5, m1, m3
vpsubd m6, m0, m2
vpsubd m7, m1, m3

vpaddd m4, m2, [esp-16]
vpaddd m5, m3, [esp]
vpsubd m6, m2, [esp-16]
vpsubd m7, m3, [esp]

vpabsd m4, m4
vpabsd m5, m5
vpabsd m6, m6
vpabsd m7, m7

vpaddd m2, m4, m5
vpaddd m2, m6
vpaddd m2, m7

vpaddd m4, m0, [esp-16*3]
vpaddd m5, m1, [esp-16*2]
vpsubd m6, m0, [esp-16*3]
vpsubd m7, m1, [esp-16*2]

vpabsd m4, m4
vpabsd m5, m5
vpabsd m6, m6
vpabsd m7, m7

vpaddd m0, m4, m5
vpaddd m0, m6
vpaddd m0, m7

vpaddd m0, m2

%endif

%if ARCH_X86_64
vpabsd m0, m0
vpabsd m1, m1
vpabsd m2, m2
vpabsd m3, m3
vpabsd m4, m4
vpabsd m5, m5
vpabsd m6, m6
vpabsd m7, m7

vpaddd m0, m1
vpaddd m0, m2
vpaddd m0, m3
vpaddd m0, m4
vpaddd m0, m5
vpaddd m0, m6
vpaddd m0, m7
%endif

vphaddd m0, m0
vphaddd m0, m0
vpextrd eax, m0, 1
vpinsrd m1, eax, 0
vpaddd m0, m1
vpextrd eax, m0, 1
vpinsrd m1, eax, 0
vpaddd m0, m1

%if ARCH_X86_64 == 0
    lea esp, [esp+16*4]
%endif

vmovd eax, m0

;Uncomment if transformed values not divided elsewhere
;add eax, 2
;shr eax, 2

RET