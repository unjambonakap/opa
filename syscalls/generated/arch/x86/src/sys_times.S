#include <opa_asm.h>
    ENTRY(times)


    movl    $__NR_times, %eax
    int     $0x80


    ret
END(times)

