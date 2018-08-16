#include <opa_asm.h>
    ENTRY(madvise)


    movl    $__NR_madvise, %eax
    syscall


    ret
END(madvise)

