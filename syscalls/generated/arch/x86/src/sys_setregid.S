#include <opa_asm.h>
    ENTRY(setregid)


    movl    $__NR_setregid32, %eax
    int     $0x80


    ret
END(setregid)

