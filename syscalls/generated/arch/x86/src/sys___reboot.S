#include <opa_asm.h>
    ENTRY(__reboot)


    movl    $__NR_reboot, %eax
    int     $0x80


    ret
END(__reboot)

