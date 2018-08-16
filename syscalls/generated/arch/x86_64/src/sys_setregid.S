#include <opa_asm.h>
    ENTRY(setregid)


    movl    $__NR_setregid, %eax
    syscall


    ret
END(setregid)

