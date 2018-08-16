#include <opa_asm.h>
    ENTRY(acct)


    movl    $__NR_acct, %eax
    int     $0x80


    ret
END(acct)

