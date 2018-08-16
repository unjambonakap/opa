#include <opa_asm.h>
    ENTRY(__llseek)


    movl    $__NR__llseek, %eax
    int     $0x80


    ret
END(__llseek)

