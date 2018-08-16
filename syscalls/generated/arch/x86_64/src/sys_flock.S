#include <opa_asm.h>
    ENTRY(flock)


    movl    $__NR_flock, %eax
    syscall


    ret
END(flock)

