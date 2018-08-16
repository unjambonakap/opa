#include <opa_asm.h>
    ENTRY(madvise)


    movl    $__NR_madvise, %eax
    int     $0x80


    ret
END(madvise)

