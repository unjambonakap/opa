#include <opa_asm.h>
    ENTRY(clock_getres)


    movl    $__NR_clock_getres, %eax
    int     $0x80


    ret
END(clock_getres)

