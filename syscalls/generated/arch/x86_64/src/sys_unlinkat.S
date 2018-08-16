#include <opa_asm.h>
    ENTRY(unlinkat)


    movl    $__NR_unlinkat, %eax
    syscall


    ret
END(unlinkat)

