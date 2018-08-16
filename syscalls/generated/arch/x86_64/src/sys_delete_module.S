#include <opa_asm.h>
    ENTRY(delete_module)


    movl    $__NR_delete_module, %eax
    syscall


    ret
END(delete_module)

