#include <opa_asm.h>
    ENTRY(init_module)


    movl    $__NR_init_module, %eax
    int     $0x80


    ret
END(init_module)

