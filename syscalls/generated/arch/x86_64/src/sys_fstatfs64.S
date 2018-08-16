#include <opa_asm.h>
    ENTRY(fstatfs64)


    movl    $__NR_fstatfs, %eax
    syscall


    ret
END(fstatfs64)

