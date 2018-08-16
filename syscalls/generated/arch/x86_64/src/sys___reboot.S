#include <opa_asm.h>
    ENTRY(__reboot)

mov %rcx, %r10

    movl    $__NR_reboot, %eax
    syscall


    ret
END(__reboot)

