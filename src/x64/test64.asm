; Function to get CPUID for identifying CPU capabilities
bits 64
section .code
global cpuId64

cpuId64:
    mov  r8, rcx ; pointer to ecx-output
    mov  r9, rdx ; pointer to edx-output

    mov eax,1
    cpuid
    mov dword [r8], ecx
    mov dword [r9], edx
    mov eax,0
    ret