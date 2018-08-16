#include <opa_asm.h>
    ENTRY(tgkill)


    movl    $__NR_tgkill, %eax
    syscall


    ret
END(tgkill)

