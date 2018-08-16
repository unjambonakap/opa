#include <opa_asm.h>
    ENTRY(swapon)


    movl    $__NR_swapon, %eax
    int     $0x80


    ret
END(swapon)

