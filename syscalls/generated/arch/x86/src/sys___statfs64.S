#include <opa_asm.h>
    ENTRY(__statfs64)


    movl    $__NR_statfs64, %eax
    int     $0x80


    ret
END(__statfs64)

