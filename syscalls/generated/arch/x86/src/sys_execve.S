#include <opa_asm.h>
    ENTRY(execve)


    movl    $__NR_execve, %eax
    int     $0x80


    ret
END(execve)

