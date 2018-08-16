#include <opa_asm.h>
    ENTRY(chroot)


    movl    $__NR_chroot, %eax
    int     $0x80


    ret
END(chroot)

