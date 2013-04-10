; Function to get CPUID for identifying CPU capabilities
bits 32
global _cpuId

_cpuId:    
    mov eax,1
    cpuid
    mov eax, dword [esp+4]
    mov dword [eax], ecx
    mov eax, dword [esp+8]
    mov dword [eax], edx
    mov eax,0
    ret