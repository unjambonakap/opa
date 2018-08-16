#include <opa_asm.h>
    ENTRY(mlockall)


    movl    $__NR_mlockall, %eax
    syscall


    ret
END(mlockall)

