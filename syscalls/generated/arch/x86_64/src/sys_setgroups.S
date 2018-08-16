#include <opa_asm.h>
    ENTRY(setgroups)


    movl    $__NR_setgroups, %eax
    syscall


    ret
END(setgroups)

