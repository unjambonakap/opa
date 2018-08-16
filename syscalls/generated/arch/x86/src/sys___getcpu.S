#include <opa_asm.h>
    ENTRY(__getcpu)


    movl    $__NR_getcpu, %eax
    int     $0x80


    ret
END(__getcpu)

