#include <opa_asm.h>
    ENTRY(symlinkat)


    movl    $__NR_symlinkat, %eax
    syscall


    ret
END(symlinkat)

