#include <opa_asm.h>
    ENTRY(dup3)


    movl    $__NR_dup3, %eax
    syscall


    ret
END(dup3)

