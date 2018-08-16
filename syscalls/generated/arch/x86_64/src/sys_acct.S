#include <opa_asm.h>
    ENTRY(acct)


    movl    $__NR_acct, %eax
    syscall


    ret
END(acct)

