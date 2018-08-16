#include <opa_asm.h>
    ENTRY(munlockall)


    movl    $__NR_munlockall, %eax
    int     $0x80


    ret
END(munlockall)

