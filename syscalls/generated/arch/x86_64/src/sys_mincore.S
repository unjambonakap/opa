#include <opa_asm.h>
    ENTRY(mincore)


    movl    $__NR_mincore, %eax
    syscall


    ret
END(mincore)

