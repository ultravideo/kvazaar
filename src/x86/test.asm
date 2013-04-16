; Function to get CPUID for identifying CPU capabilities
bits 32
section .code
global cpuId32

;void __cdecl cpuId(int* ecx, int *edx );

cpuId32:    
    mov eax,1
    cpuid
    mov eax, dword [esp+4] ; pointer to ecx-output, pushed second to the stack
    mov dword [eax], ecx
    mov eax, dword [esp+8] ; pointer to edx-output, pushed first to the stack
    mov dword [eax], edx
    mov eax,0
    ret