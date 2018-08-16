#include <opa_asm.h>
    ENTRY(mlockall)


    movl    $__NR_mlockall, %eax
    int     $0x80


    ret
END(mlockall)

