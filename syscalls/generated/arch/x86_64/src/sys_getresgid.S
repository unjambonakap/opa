#include <opa_asm.h>
    ENTRY(getresgid)


    movl    $__NR_getresgid, %eax
    syscall


    ret
END(getresgid)

