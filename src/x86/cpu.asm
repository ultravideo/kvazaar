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

; Function to get CPUID for identifying CPU capabilities

;void cpuid(int *ecx, int *edx);
SECTION .text
; we allocate 5 registers just because R2-R3 can be registers we need from cpuid (ecx,edx)
cglobal cpu_cpuid, 2,5
  ;rbx is not caller saved and cpuid modifies it
  push rbx
  ;Depending on the architecture cpuid might modify r0 and r1 so store them in stack
  push r0
  push r1
  mov eax,1
  cpuid
  ; pop out r1, the second parameter,  containing pointer to "edx" data
  pop r4
  mov dword [r4], edx ; store 32bit edx value to the pointer
  ; pop out r0, the first parameter, containing pointer to "ecx" data
  pop r4
  mov dword [r4], ecx ; store 32bit ecx value to the pointer
  pop rbx
  ; library handles pushing and poping of the needed registers with RET macro
  RET