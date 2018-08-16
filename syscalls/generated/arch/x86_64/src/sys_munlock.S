#include <opa_asm.h>
    ENTRY(munlock)


    movl    $__NR_munlock, %eax
    syscall


    ret
END(munlock)

