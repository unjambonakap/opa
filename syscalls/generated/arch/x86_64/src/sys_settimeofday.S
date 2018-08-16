#include <opa_asm.h>
    ENTRY(settimeofday)


    movl    $__NR_settimeofday, %eax
    syscall


    ret
END(settimeofday)

