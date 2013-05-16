; Function to get CPUID for identifying CPU capabilities
bits 32
section .code
global _cpuId32

;void __cdecl cpuId(int* ecx, int *edx );

_cpuId32:    
    push eax
    push ecx
    push edx
    mov eax,1
    cpuid
    mov eax, dword [esp+4] ; pointer to ecx-output, pushed second to the stack
    mov dword [eax], ecx
    mov eax, dword [esp+8] ; pointer to edx-output, pushed first to the stack
    mov dword [eax], edx
    pop edx
    pop ecx
    pop eax
    ret