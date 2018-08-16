#include <opa_asm.h>
    ENTRY(__epoll_pwait)

mov %rcx, %r10

    movl    $__NR_epoll_pwait, %eax
    syscall


    ret
END(__epoll_pwait)

