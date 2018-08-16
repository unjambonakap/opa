#include <opa_asm.h>
    ENTRY(recvmsg)


    movl    $__NR_recvmsg, %eax
    syscall


    ret
END(recvmsg)

