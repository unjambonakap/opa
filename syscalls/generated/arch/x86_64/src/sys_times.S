#include <opa_asm.h>
    ENTRY(times)


    movl    $__NR_times, %eax
    syscall


    ret
END(times)

