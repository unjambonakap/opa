#include <opa_asm.h>
    ENTRY(__pselect6)


    movl    $__NR_pselect6, %eax
    int     $0x80


    ret
END(__pselect6)

