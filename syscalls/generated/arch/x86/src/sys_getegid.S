#include <opa_asm.h>
    ENTRY(getegid)


    movl    $__NR_getegid32, %eax
    int     $0x80


    ret
END(getegid)

