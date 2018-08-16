#include <opa_asm.h>
    ENTRY(personality)


    movl    $__NR_personality, %eax
    syscall


    ret
END(personality)

