#include <opa_asm.h>
    ENTRY(uname)


    movl    $__NR_uname, %eax
    int     $0x80


    ret
END(uname)

