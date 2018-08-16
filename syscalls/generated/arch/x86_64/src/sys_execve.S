#include <opa_asm.h>
    ENTRY(execve)


    movl    $__NR_execve, %eax
    syscall


    ret
END(execve)

