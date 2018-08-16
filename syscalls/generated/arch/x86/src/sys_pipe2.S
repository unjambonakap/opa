#include <opa_asm.h>
    ENTRY(pipe2)


    movl    $__NR_pipe2, %eax
    int     $0x80


    ret
END(pipe2)

