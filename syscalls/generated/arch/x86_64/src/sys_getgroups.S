#include <opa_asm.h>
    ENTRY(getgroups)


    movl    $__NR_getgroups, %eax
    syscall


    ret
END(getgroups)

