#include <opa_asm.h>
    ENTRY(__gettimeofday)


    movl    $__NR_gettimeofday, %eax
    syscall


    ret
END(__gettimeofday)

