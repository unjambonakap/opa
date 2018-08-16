#include <opa_asm.h>
    ENTRY(mlock)


    movl    $__NR_mlock, %eax
    int     $0x80


    ret
END(mlock)

