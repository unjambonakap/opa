#include <opa_asm.h>
    ENTRY(__clock_nanosleep)


    movl    $__NR_clock_nanosleep, %eax
    int     $0x80


    ret
END(__clock_nanosleep)

