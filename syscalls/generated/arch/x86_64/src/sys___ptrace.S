#include <opa_asm.h>
    ENTRY(__ptrace)

mov %rcx, %r10

    movl    $__NR_ptrace, %eax
    syscall


    ret
END(__ptrace)

