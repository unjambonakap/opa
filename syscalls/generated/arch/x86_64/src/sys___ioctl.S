#include <opa_asm.h>
    ENTRY(__ioctl)


    movl    $__NR_ioctl, %eax
    syscall


    ret
END(__ioctl)

