#include <opa_asm.h>
    ENTRY(unshare)


    movl    $__NR_unshare, %eax
    int     $0x80


    ret
END(unshare)

