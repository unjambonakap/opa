#include <opa_asm.h>
    ENTRY(getgid)


    movl    $__NR_getgid32, %eax
    int     $0x80


    ret
END(getgid)

