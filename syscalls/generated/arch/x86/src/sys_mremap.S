#include <opa_asm.h>
    ENTRY(mremap)


    movl    $__NR_mremap, %eax
    int     $0x80


    ret
END(mremap)

