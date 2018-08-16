#include <opa_asm.h>
    ENTRY(umask)


    movl    $__NR_umask, %eax
    syscall


    ret
END(umask)

