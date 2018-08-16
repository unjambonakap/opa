#include <opa_asm.h>
    ENTRY(__fstatfs64)


    movl    $__NR_fstatfs64, %eax
    int     $0x80


    ret
END(__fstatfs64)

