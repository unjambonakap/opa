#include <opa_asm.h>
    ENTRY(setresgid)


    movl    $__NR_setresgid, %eax
    syscall


    ret
END(setresgid)

