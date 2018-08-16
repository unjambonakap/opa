#include <opa_asm.h>
    ENTRY(fgetxattr)


    movl    $__NR_fgetxattr, %eax
    int     $0x80


    ret
END(fgetxattr)

