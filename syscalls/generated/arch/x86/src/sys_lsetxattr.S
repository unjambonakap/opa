#include <opa_asm.h>
    ENTRY(lsetxattr)


    movl    $__NR_lsetxattr, %eax
    int     $0x80


    ret
END(lsetxattr)

