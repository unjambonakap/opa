#include <opa_asm.h>
    ENTRY(init_module)


    movl    $__NR_init_module, %eax
    syscall


    ret
END(init_module)

