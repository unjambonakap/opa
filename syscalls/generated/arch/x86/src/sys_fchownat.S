#include <opa_asm.h>
    ENTRY(fchownat)


    movl    $__NR_fchownat, %eax
    int     $0x80


    ret
END(fchownat)

