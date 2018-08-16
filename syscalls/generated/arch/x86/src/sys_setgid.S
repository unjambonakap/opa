#include <opa_asm.h>
    ENTRY(setgid)


    movl    $__NR_setgid32, %eax
    int     $0x80


    ret
END(setgid)

