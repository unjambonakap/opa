#include <opa_asm.h>
    ENTRY(clock_getres)


    movl    $__NR_clock_getres, %eax
    syscall


    ret
END(clock_getres)

