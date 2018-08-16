#include <opa_asm.h>
    ENTRY(read)


    movl    $__NR_read, %eax
    int     $0x80


    ret
END(read)

