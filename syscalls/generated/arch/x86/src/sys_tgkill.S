#include <opa_asm.h>
    ENTRY(tgkill)


    movl    $__NR_tgkill, %eax
    int     $0x80


    ret
END(tgkill)

