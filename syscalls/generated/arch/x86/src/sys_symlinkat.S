#include <opa_asm.h>
    ENTRY(symlinkat)


    movl    $__NR_symlinkat, %eax
    int     $0x80


    ret
END(symlinkat)

