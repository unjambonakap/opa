#include <opa_asm.h>
    ENTRY(umount2)


    movl    $__NR_umount2, %eax
    int     $0x80


    ret
END(umount2)

