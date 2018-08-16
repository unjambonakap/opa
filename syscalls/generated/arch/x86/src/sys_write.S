#include <opa_asm.h>
    ENTRY(write)


    movl    $__NR_write, %eax
    int     $0x80


    ret
END(write)

