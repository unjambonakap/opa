#include <opa_asm.h>
    ENTRY(setpgid)


    movl    $__NR_setpgid, %eax
    int     $0x80


    ret
END(setpgid)

