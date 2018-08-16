#include <opa_asm.h>
    ENTRY(swapon)


    movl    $__NR_swapon, %eax
    syscall


    ret
END(swapon)

