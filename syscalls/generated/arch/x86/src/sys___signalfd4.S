#include <opa_asm.h>
    ENTRY(__signalfd4)


    movl    $__NR_signalfd4, %eax
    int     $0x80


    ret
END(__signalfd4)

