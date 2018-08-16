#include <opa_asm.h>
    ENTRY(pread64)


    movl    $__NR_pread64, %eax
    int     $0x80


    ret
END(pread64)

