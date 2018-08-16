#include <opa_asm.h>
    ENTRY(fsetxattr)


    movl    $__NR_fsetxattr, %eax
    int     $0x80


    ret
END(fsetxattr)

