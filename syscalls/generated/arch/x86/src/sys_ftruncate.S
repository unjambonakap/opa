#include <opa_asm.h>
    ENTRY(ftruncate)


    movl    $__NR_ftruncate, %eax
    int     $0x80


    ret
END(ftruncate)

