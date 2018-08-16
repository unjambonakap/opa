#include <opa_asm.h>
    ENTRY(chroot)


    movl    $__NR_chroot, %eax
    syscall


    ret
END(chroot)

