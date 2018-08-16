#include <opa_asm.h>
    ENTRY(clock_settime)


    movl    $__NR_clock_settime, %eax
    syscall


    ret
END(clock_settime)

