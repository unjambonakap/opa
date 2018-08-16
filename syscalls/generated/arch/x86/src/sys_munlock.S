#include <opa_asm.h>
    ENTRY(munlock)


    movl    $__NR_munlock, %eax
    int     $0x80


    ret
END(munlock)

