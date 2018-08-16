#include <opa_asm.h>
    ENTRY(fstat64)


    movl    $__NR_fstat, %eax
    syscall


    ret
END(fstat64)

