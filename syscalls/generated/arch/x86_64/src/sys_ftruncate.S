#include <opa_asm.h>
    ENTRY(ftruncate)


    movl    $__NR_ftruncate, %eax
    syscall


    ret
END(ftruncate)

