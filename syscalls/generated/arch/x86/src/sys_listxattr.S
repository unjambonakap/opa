#include <opa_asm.h>
    ENTRY(listxattr)


    movl    $__NR_listxattr, %eax
    int     $0x80


    ret
END(listxattr)

