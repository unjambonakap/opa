#include <opa_asm.h>
    ENTRY(mknodat)


    movl    $__NR_mknodat, %eax
    int     $0x80


    ret
END(mknodat)

