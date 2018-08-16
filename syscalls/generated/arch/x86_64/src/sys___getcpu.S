#include <opa_asm.h>
    ENTRY(__getcpu)


    movl    $__NR_getcpu, %eax
    syscall


    ret
END(__getcpu)

