#include <opa_asm.h>
    ENTRY(mount)


    movl    $__NR_mount, %eax
    int     $0x80


    ret
END(mount)

