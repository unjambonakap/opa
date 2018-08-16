#include <opa_asm.h>
    ENTRY(capset)


    movl    $__NR_capset, %eax
    syscall


    ret
END(capset)

