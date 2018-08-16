#include <opa_asm.h>
    ENTRY(sendmsg)


    movl    $__NR_sendmsg, %eax
    syscall


    ret
END(sendmsg)

