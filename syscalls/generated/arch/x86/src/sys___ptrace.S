#include <opa_asm.h>
    ENTRY(__ptrace)


    movl    $__NR_ptrace, %eax
    int     $0x80


    ret
END(__ptrace)

