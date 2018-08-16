#include <opa_asm.h>
    ENTRY(uname)


    movl    $__NR_uname, %eax
    syscall


    ret
END(uname)

