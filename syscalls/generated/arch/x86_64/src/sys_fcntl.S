#include <opa_asm.h>
    ENTRY(fcntl)


    movl    $__NR_fcntl, %eax
    syscall


    ret
END(fcntl)

