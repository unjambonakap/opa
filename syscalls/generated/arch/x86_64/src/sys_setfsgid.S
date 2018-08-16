#include <opa_asm.h>
    ENTRY(setfsgid)


    movl    $__NR_setfsgid, %eax
    syscall


    ret
END(setfsgid)

