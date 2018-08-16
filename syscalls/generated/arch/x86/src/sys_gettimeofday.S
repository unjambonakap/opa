#include <opa_asm.h>
    ENTRY(gettimeofday)


    movl    $__NR_gettimeofday, %eax
    int     $0x80


    ret
END(gettimeofday)

