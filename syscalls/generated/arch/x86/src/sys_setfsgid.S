#include <opa_asm.h>
    ENTRY(setfsgid)


    movl    $__NR_setfsgid, %eax
    int     $0x80


    ret
END(setfsgid)

