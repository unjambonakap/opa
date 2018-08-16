#include <opa_asm.h>
    ENTRY(__clock_gettime)


    movl    $__NR_clock_gettime, %eax
    syscall


    ret
END(__clock_gettime)

