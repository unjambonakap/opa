#include <opa_asm.h>
    ENTRY(read)


    movl    $__NR_read, %eax
    syscall


    ret
END(read)

