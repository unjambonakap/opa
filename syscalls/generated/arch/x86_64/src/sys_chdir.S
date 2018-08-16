#include <opa_asm.h>
    ENTRY(chdir)


    movl    $__NR_chdir, %eax
    syscall


    ret
END(chdir)

