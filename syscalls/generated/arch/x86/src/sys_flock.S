#include <opa_asm.h>
    ENTRY(flock)


    movl    $__NR_flock, %eax
    int     $0x80


    ret
END(flock)

