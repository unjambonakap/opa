#include <opa_asm.h>
    ENTRY(settimeofday)


    movl    $__NR_settimeofday, %eax
    int     $0x80


    ret
END(settimeofday)

