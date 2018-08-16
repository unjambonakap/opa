#include <opa_asm.h>
    ENTRY(__fadvise64)


    movl    $__NR_fadvise64_64, %eax
    int     $0x80


    ret
END(__fadvise64)

