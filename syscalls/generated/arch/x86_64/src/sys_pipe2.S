#include <opa_asm.h>
    ENTRY(pipe2)


    movl    $__NR_pipe2, %eax
    syscall


    ret
END(pipe2)

