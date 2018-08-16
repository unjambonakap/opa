#include <opa_asm.h>
    ENTRY(getresgid)


    movl    $__NR_getresgid32, %eax
    int     $0x80


    ret
END(getresgid)

