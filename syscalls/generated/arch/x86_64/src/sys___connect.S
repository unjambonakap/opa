#include <opa_asm.h>
    ENTRY(__connect)


    movl    $__NR_connect, %eax
    syscall


    ret
END(__connect)

