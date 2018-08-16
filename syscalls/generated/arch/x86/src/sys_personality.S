#include <opa_asm.h>
    ENTRY(personality)


    movl    $__NR_personality, %eax
    int     $0x80


    ret
END(personality)

