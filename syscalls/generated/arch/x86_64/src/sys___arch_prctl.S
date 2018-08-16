#include <opa_asm.h>
    ENTRY(__arch_prctl)


    movl    $__NR_arch_prctl, %eax
    syscall


    ret
END(__arch_prctl)

