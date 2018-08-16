#include <opa_asm.h>
    ENTRY(getppid)


    movl    $__NR_getppid, %eax
    int     $0x80


    ret
END(getppid)

