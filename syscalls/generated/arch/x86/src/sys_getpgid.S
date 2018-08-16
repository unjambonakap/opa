#include <opa_asm.h>
    ENTRY(getpgid)


    movl    $__NR_getpgid, %eax
    int     $0x80


    ret
END(getpgid)

