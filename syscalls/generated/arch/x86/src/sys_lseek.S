#include <opa_asm.h>
    ENTRY(lseek)


    movl    $__NR_lseek, %eax
    int     $0x80


    ret
END(lseek)

