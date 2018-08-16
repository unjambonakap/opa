#include <opa_asm.h>
    ENTRY(lseek)


    movl    $__NR_lseek, %eax
    syscall


    ret
END(lseek)

