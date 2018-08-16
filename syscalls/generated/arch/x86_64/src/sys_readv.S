#include <opa_asm.h>
    ENTRY(readv)


    movl    $__NR_readv, %eax
    syscall


    ret
END(readv)

