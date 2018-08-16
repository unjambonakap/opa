#include <opa_asm.h>
    ENTRY(__waitid)


    movl    $__NR_waitid, %eax
    int     $0x80


    ret
END(__waitid)

