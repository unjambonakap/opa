#include <opa_asm.h>
    ENTRY(clock_gettime)


    movl    $__NR_clock_gettime, %eax
    int     $0x80


    ret
END(clock_gettime)

