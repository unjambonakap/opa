#include <opa_asm.h>
    ENTRY(prlimit64)


    movl    $__NR_prlimit64, %eax
    int     $0x80


    ret
END(prlimit64)

