#include <opa_asm.h>
    ENTRY(close)


    movl    $__NR_close, %eax
    syscall


    ret
END(close)

