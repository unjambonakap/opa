#include <opa_asm.h>
    ENTRY(dup)


    movl    $__NR_dup, %eax
    syscall


    ret
END(dup)

