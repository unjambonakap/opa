#include <opa_asm.h>
    ENTRY(sethostname)


    movl    $__NR_sethostname, %eax
    syscall


    ret
END(sethostname)

