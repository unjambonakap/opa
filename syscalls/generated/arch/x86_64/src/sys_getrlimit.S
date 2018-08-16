#include <opa_asm.h>
    ENTRY(getrlimit)


    movl    $__NR_getrlimit, %eax
    syscall


    ret
END(getrlimit)

