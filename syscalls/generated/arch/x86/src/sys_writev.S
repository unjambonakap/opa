#include <opa_asm.h>
    ENTRY(writev)


    movl    $__NR_writev, %eax
    int     $0x80


    ret
END(writev)

