#include <opa_asm.h>
    ENTRY(shutdown)


    movl    $__NR_shutdown, %eax
    syscall


    ret
END(shutdown)

