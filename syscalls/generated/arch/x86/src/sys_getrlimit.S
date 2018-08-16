#include <opa_asm.h>
    ENTRY(getrlimit)


    movl    $__NR_ugetrlimit, %eax
    int     $0x80


    ret
END(getrlimit)

