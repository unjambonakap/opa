#include <opa_asm.h>
    ENTRY(getgroups)


    movl    $__NR_getgroups32, %eax
    int     $0x80


    ret
END(getgroups)

