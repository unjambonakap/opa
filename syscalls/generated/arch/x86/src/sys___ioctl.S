#include <opa_asm.h>
    ENTRY(__ioctl)


    movl    $__NR_ioctl, %eax
    int     $0x80


    ret
END(__ioctl)

