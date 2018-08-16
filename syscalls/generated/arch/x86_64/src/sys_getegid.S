#include <opa_asm.h>
    ENTRY(getegid)


    movl    $__NR_getegid, %eax
    syscall


    ret
END(getegid)

