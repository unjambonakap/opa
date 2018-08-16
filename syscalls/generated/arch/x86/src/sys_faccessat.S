#include <opa_asm.h>
    ENTRY(faccessat)


    movl    $__NR_faccessat, %eax
    int     $0x80


    ret
END(faccessat)

