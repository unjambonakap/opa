#include <opa_asm.h>
    ENTRY(clock_settime)


    movl    $__NR_clock_settime, %eax
    int     $0x80


    ret
END(clock_settime)

