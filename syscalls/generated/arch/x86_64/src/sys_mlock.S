#include <opa_asm.h>
    ENTRY(mlock)


    movl    $__NR_mlock, %eax
    syscall


    ret
END(mlock)

