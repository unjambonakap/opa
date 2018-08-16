#include <opa_asm.h>
    ENTRY(__ppoll)


    movl    $__NR_ppoll, %eax
    int     $0x80


    ret
END(__ppoll)

