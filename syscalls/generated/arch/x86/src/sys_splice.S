#include <opa_asm.h>
    ENTRY(splice)


    movl    $__NR_splice, %eax
    int     $0x80


    ret
END(splice)

