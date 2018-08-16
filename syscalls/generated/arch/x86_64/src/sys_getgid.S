#include <opa_asm.h>
    ENTRY(getgid)


    movl    $__NR_getgid, %eax
    syscall


    ret
END(getgid)

