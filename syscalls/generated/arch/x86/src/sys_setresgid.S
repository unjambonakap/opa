#include <opa_asm.h>
    ENTRY(setresgid)


    movl    $__NR_setresgid32, %eax
    int     $0x80


    ret
END(setresgid)

