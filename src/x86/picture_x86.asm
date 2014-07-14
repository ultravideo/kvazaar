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

;cglobal and RET macros are from the x86.inc
;they push and pop the necessary registers to
;stack depending on the operating system

;Usage: cglobal name, %1, %2, %3
;1%: Number of arguments
;2%: Number of registers used
;3%: Number of xmm registers used.
;More info in x86inc.asm

SECTION .text

;Set x86inc.asm macros to use avx and xmm registers
INIT_XMM avx

;KVZ_SAD_4X4
;Calculates SAD of the 16 consequtive bytes in memory
;r0 address of the first value(current frame)
;r1 address of the first value(reference frame)

cglobal sad_4x4, 2, 2, 2

;Load 16 bytes of both frames
vmovdqu m0, [r0]
vmovdqu m1, [r1]

;Calculate SAD. The results are written in
;m0[15:0] and m0[79:64]
vpsadbw m0, m1

;Sum the result
vmovhlps m1, m0
vpaddw m0, m1

;Set the result to eax
vmovd eax, m0

RET


;KVZ_SAD_4X4_STRIDE
;Calculates SAD of a 4x4 block inside a frame with stride
;r0 address of the first value(current)
;r1 address of the first value(reference)
;r2 stride

cglobal sad_4x4_stride, 3, 3, 2

;Load 4 times 4 bytes of both frames
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

;Load the first half of both frames
vmovdqu m0, [r0]
vmovdqu m2, [r0+16]

vmovdqu m1, [r1]
vmovdqu m3, [r1+16]

;Calculate SADs for both
vpsadbw m0, m1
vpsadbw m2, m3

;Sum
vpaddw m0, m2

;Repeat for the latter half
vmovdqu m1, [r0+16*2]
vmovdqu m3, [r0+16*3]

vmovdqu m2, [r1+16*2]
vmovdqu m4, [r1+16*3]

vpsadbw m1, m2
vpsadbw m3, m4

vpaddw m1, m3

;Sum all the SADs
vpaddw m0, m1

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

;Zero m0 register
vpxor m0, m0

;Load the first half to m1 and m3 registers(cur)
;Current frame
;Load to the high 64 bits of xmm
vmovhpd m1, [r0]
add r0, r2
;Load to the low 64 bits
vmovlpd m1, [r0] 

vmovhpd m3, [r0+r2]
vmovlpd m3, [r0+r2*2] 
;lea calculates the address to r0,
;but doesn't load anything from
;the memory. Equivalent for
;two add r0, r2 instructions.
lea r0, [r0+r2*2]
add r0, r2

;Reference frame
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

;Repeat for the other half
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

;Zero m4
vpxor m4, m4

%assign i 0

;Repeat 8 times.
%rep 8

;Load the next to rows of the current frame
vmovdqu m0, [r0 + 16 * i]
vmovdqu m2, [r0 + 16 * (i + 1)]

;Load the next to rows of the reference frame
vmovdqu m1, [r1 + 16 * i]
vmovdqu m3, [r1 + 16 * (i + 1)]

vpsadbw m0, m1
vpsadbw m2, m3

;Accumulate SADs to m4
vpaddw m4, m0
vpaddw m4, m2

%assign i i+2

%endrep

;Calculate the final sum
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

; Load the next 2 rows from rec_buf to m0 and m2
vmovdqu m0, [r0]
vmovdqu m2, [r0 + r2]
lea r0, [r0 + r2*2]

; Load the next 2 rows from ref_buf to m1 and m3
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

;Load 8 bytes from memory and zero extend
;to 16-bit values. Calculate difference.
vpmovzxbw m0, [r0]
vpmovzxbw m2, [r1]
vpsubw m0, m2

vpmovzxbw m1, [r0+8]
vpmovzxbw m3, [r1+8]
vpsubw m1, m3

;Hadamard transform
;Horizontal phase
vphaddw m4, m0, m1
vphsubw m5, m0, m1

vphaddw m0, m4, m5
vphsubw m1, m4, m5

;Vertical phase
vphaddw m4, m0, m1
vphsubw m5, m0, m1

vphaddw m0, m4, m5
vphsubw m1, m4, m5

;Calculate absolute values
vpabsw m0, m0
vpabsw m1, m1

;Sum the all the transformed values
vpaddw m0, m1

vphaddw m0, m0
vphaddw m0, m0
vphaddw m0, m0

;Extract the lowest 16 bits of m0
;into eax
vpextrw eax, m0, 0

;Add 1 and divide by 2
add eax, 1
shr eax, 1

RET

;KVZ_ZERO_EXTEND
;zero extend all packed words in xmm to dwords in ymm
;%1 number of the destination register
;%2 number of a free xmm register
%macro KVZ_ZERO_EXTEND 2

;Zero extend low 64 bits
vpmovzxwd xmm%2, xmm%1
;Zero extend high 64 bits
vmovhlps xmm%1, xmm%1
vpmovzxwd xmm%1, xmm%1
;Combine xmm%1 and xmm%2 in ymm%1
vinserti128 ymm%1, ymm%2, xmm%1, 1

%endmacro

;KVZ_SATD_8X8_STRIDE
;Calculates SATD of a 8x8 block inside a frame with stride
;r0 address of the first value(reference)
;r1 address of the first value(current)
;r2 stride
;
;The Result is written in the register r4

%macro KVZ_SATD_8X8_STRIDE 0


;Calculate differences of the 8 rows into
;registers m0-m7
vpmovzxbw m0, [r0]
vpmovzxbw m7, [r2]
vpsubw m0, m7

vpmovzxbw m1, [r0+r1]
vpmovzxbw m7, [r2+r3]
vpsubw m1, m7

;Set r0 and r2 2 rows forward
lea r0, [r0+r1*2]
lea r2, [r2+r3*2]

vpmovzxbw m2, [r0]
vpmovzxbw m7, [r2]
vpsubw m2, m7

vpmovzxbw m3, [r0+r1]
vpmovzxbw m7, [r2+r3]
vpsubw m3, m7

lea r0, [r0+r1*2]
lea r2, [r2+r3*2]

vpmovzxbw m4, [r0]
vpmovzxbw m7, [r2]
vpsubw m4, m7

vpmovzxbw m5, [r0+r1]
vpmovzxbw m7, [r2+r3]
vpsubw m5, m7

lea r0, [r0+r1*2]
lea r2, [r2+r3*2]

vpmovzxbw m6, [r0]
vpmovzxbw m7, [r2]
vpsubw m6, m7

;32-bit AVX doesn't have registers
;xmm8-xmm15, use stack instead
%if ARCH_X86_64
    vpmovzxbw m7, [r0+r1]
    vpmovzxbw m8, [r2+r3]
    vpsubw m7, m8
%else
    %define temp0 esp+16*3
    %define temp1 esp+16*2
    %define temp2 esp+16*1
    %define temp3 esp+16*0
    ;Reserve memory for 4 x 128 bits
    sub esp, 16*4

    vpmovzxbw m7, [r2+r3]
    vmovdqu [temp0], m7
    vpmovzxbw m7, [r0+r1]
    vpsubw m7, [temp0]

    ;Put rows 5-8 to stack
    vmovdqu [temp0], m4
    vmovdqu [temp1], m5
    vmovdqu [temp2], m6
    vmovdqu [temp3], m7
%endif

;Hadamard transform
;Horizontal phaze

%if ARCH_X86_64
;Calculate each step into the other half of the
;xmm registers(xmm0-xmm7 <-> xmm8-xmm15) in order to
;eliminate the need for memory access.
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

%else
;Calculate horizontal transforms for the first half.
;Then load to second half into registers and store
;transforms in the stack.
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


    vmovdqu m3, [temp3]
    vmovdqu m2, [temp2]
    vmovdqu m1, [temp1]
    vmovdqu m0, [temp0]

    vmovdqu [temp3], m7
    vmovdqu [temp2], m6
    vmovdqu [temp1], m5
    vmovdqu [temp0], m4

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

    vpaddw m8, m0, m2
    vpaddw m9, m1, m3
    vpsubw m10, m0, m2
    vpsubw m11, m1, m3

    vpaddw m12, m4, m6
    vpaddw m13, m5, m7
    vpsubw m14, m4, m6
    vpsubw m15, m5, m7

    vpaddw m0, m8, m12
    vpaddw m1, m9, m13
    vpaddw m2, m10, m14
    vpaddw m3, m11, m15

    vpsubw m4, m8, m12
    vpsubw m5, m9, m13
    vpsubw m6, m10, m14
    vpsubw m7, m11, m15

%else

    vphaddw m0, m4, m5
    vphsubw m1, m4, m5

    vphaddw m2, m6, m7
    vphsubw m3, m6, m7

    vpaddw m4, m0, m2
    vpaddw m5, m1, m3
    vpsubw m6, m0, m2
    vpsubw m7, m1, m3

    vmovdqu m3, [temp3]
    vmovdqu m2, [temp2]
    vmovdqu m1, [temp1]
    vmovdqu m0, [temp0]

    vmovdqu [temp3], m7
    vmovdqu [temp2], m6
    vmovdqu [temp1], m5
    vmovdqu [temp0], m4

    vphaddw m4, m0, m1
    vphsubw m5, m0, m1

    vphaddw m6, m2, m3
    vphsubw m7, m2, m3

    vpaddw m0, m4, m6
    vpaddw m1, m5, m7
    vpsubw m2, m4, m6
    vpsubw m3, m5, m7

    vpaddw m4, m0, [temp0]
    vpaddw m5, m1, [temp1]
    vpsubw m6, m0, [temp0]
    vpsubw m7, m1, [temp1]

    ;Calculate the absolute values and
    ;zero extend 16-bit values to 32-bit
    ;values. Needs ymm registers (256 bits).
    vpabsw m4, m4
    KVZ_ZERO_EXTEND 4, 1
    vpabsw m5, m5
    KVZ_ZERO_EXTEND 5, 1
    vpabsw m6, m6
    KVZ_ZERO_EXTEND 6, 1
    vpabsw m7, m7
    KVZ_ZERO_EXTEND 7, 1

    ;Sum half of the packed results to ymm0
    vpaddd ymm0, ymm4, ymm5
    vpaddd ymm0, ymm6
    vpaddd ymm0, ymm7

    vpaddw m4, m2, [temp2]
    vpaddw m5, m3, [temp3]
    vpsubw m6, m2, [temp2]
    vpsubw m7, m3, [temp3]

    vpabsw m4, m4
    KVZ_ZERO_EXTEND 4, 1
    vpabsw m5, m5
    KVZ_ZERO_EXTEND 5, 1
    vpabsw m6, m6
    KVZ_ZERO_EXTEND 6, 1
    vpabsw m7, m7
    KVZ_ZERO_EXTEND 7, 1

    ;Sum the rest of the packed results in ymm0
    vpaddd ymm4, ymm5
    vpaddd ymm4, ymm6
    vpaddd ymm4, ymm7

    vpaddd ymm0, ymm4

%endif

%if ARCH_X86_64
    vpabsw m0, m0
    KVZ_ZERO_EXTEND 0, 8
    vpabsw m1, m1
    KVZ_ZERO_EXTEND 1, 8
    vpabsw m2, m2
    KVZ_ZERO_EXTEND 2, 8
    vpabsw m3, m3
    KVZ_ZERO_EXTEND 3, 8
    vpabsw m4, m4
    KVZ_ZERO_EXTEND 4, 8
    vpabsw m5, m5
    KVZ_ZERO_EXTEND 5, 8
    vpabsw m6, m6
    KVZ_ZERO_EXTEND 6, 8
    vpabsw m7, m7
    KVZ_ZERO_EXTEND 7, 8

    vpaddd ymm0, ymm1
    vpaddd ymm0, ymm2
    vpaddd ymm0, ymm3
    vpaddd ymm0, ymm4
    vpaddd ymm0, ymm5
    vpaddd ymm0, ymm6
    vpaddd ymm0, ymm7
%endif

;Sum the packed values
vextracti128 m1, ymm0, 1
vpaddd m0, m1
vphaddd m0, m0
vphaddd m0, m0

;The result is in the lowest 32 bits in m0
vmovd r4d, m0

;8x8 Hadamard transform requires
;adding 2 and dividing by 4
add r4, 2
shr r4, 2

;Zero high 128 bits of ymm registers to
;prevent AVX-SSE transition penalty.
vzeroupper

%if ARCH_X86_64 == 0
    add esp, 16*4
%endif

%endmacro

;KVZ_SATD_8X8
;Calculates SATD of a 8x8 block inside a frame with stride
;r0 address of the first value(reference)
;r1 address of the first value(current)
;r2 stride
%if ARCH_X86_64
    cglobal satd_8x8, 4, 5, 16
%else
    cglobal satd_8x8, 4, 5, 8
%endif
;Set arguments
mov r2, r1
mov r1, 8
mov r3, 8
;Calculate 8x8 SATD. Result is written
;in the register r4.
KVZ_SATD_8X8_STRIDE
mov rax, r4
RET

;KVZ_SATD_NXN
;Calculates SATD of a NxN block inside a frame with stride
;r0 address of the first value(reference)
;r1 address of the first value(current)

%macro KVZ_SATD_NXN 1

%if ARCH_X86_64
    cglobal satd_%1x%1, 2, 7, 16
%else
    cglobal satd_%1x%1, 2, 7, 8
%endif
;Set arguments
mov r2, r1
mov r1, %1
mov r3, %1
;Zero r5 and r6
xor r5, r5
xor r6, r6

;Calculate SATDs of each 8x8 sub-blocks
;and accumulate the results in r6. Repeat yloop
;N times. Repeat xloop N times. r4 and r5 are counters
;for the loops.
.yloop
    ;zero r4
    xor r4, r4

    .xloop
        push r4
        
        ;Calculate SATD of the sub-block. Result is
        ;written in the register r4.
        KVZ_SATD_8X8_STRIDE
        add r6, r4

        ;Set r2 and r0 to the next sub-block
        ;on the same row
        sub r2, 6*%1-8
        sub r0, 6*%1-8

        pop r4
        add r4, 8
        cmp r4, %1
    jne .xloop

    ;Set r2 and r0 to the first sub-block
    ;on the next row(of 8x8 sub-blocks)
    add r2, 7*%1
    add r0, 7*%1

    add r5, 8
    cmp r5, %1
jne .yloop


mov rax, r6

RET

%endmacro

KVZ_SATD_NXN 16
KVZ_SATD_NXN 32
KVZ_SATD_NXN 64