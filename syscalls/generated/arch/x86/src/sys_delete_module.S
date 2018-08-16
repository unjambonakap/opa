#include <opa_asm.h>
    ENTRY(delete_module)


    movl    $__NR_delete_module, %eax
    int     $0x80


    ret
END(delete_module)

