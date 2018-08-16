#include <opa_asm.h>
    ENTRY(__brk)


    movl    $__NR_brk, %eax
    syscall


    ret
END(__brk)

