#include <opa_asm.h>
    ENTRY(fchdir)


    movl    $__NR_fchdir, %eax
    syscall


    ret
END(fchdir)

