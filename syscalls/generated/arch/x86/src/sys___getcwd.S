#include <opa_asm.h>
    ENTRY(__getcwd)


    movl    $__NR_getcwd, %eax
    int     $0x80


    ret
END(__getcwd)

