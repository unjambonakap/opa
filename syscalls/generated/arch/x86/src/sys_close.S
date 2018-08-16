#include <opa_asm.h>
    ENTRY(close)


    movl    $__NR_close, %eax
    int     $0x80


    ret
END(close)

