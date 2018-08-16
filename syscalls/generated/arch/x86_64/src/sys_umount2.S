#include <opa_asm.h>
    ENTRY(umount2)


    movl    $__NR_umount2, %eax
    syscall


    ret
END(umount2)

