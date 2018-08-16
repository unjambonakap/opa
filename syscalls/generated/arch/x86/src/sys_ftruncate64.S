#include <opa_asm.h>
    ENTRY(ftruncate64)


    movl    $__NR_ftruncate64, %eax
    int     $0x80


    ret
END(ftruncate64)

