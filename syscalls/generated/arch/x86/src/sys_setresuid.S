#include <opa_asm.h>
    ENTRY(setresuid)


    movl    $__NR_setresuid32, %eax
    int     $0x80


    ret
END(setresuid)

