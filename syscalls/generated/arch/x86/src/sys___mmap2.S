#include <opa_asm.h>
    ENTRY(__mmap2)


    movl    $__NR_mmap2, %eax
    int     $0x80


    ret
END(__mmap2)

