#include <opa_asm.h>
    ENTRY(getpgid)


    movl    $__NR_getpgid, %eax
    syscall


    ret
END(getpgid)

