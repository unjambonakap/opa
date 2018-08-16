#include <opa_asm.h>
    ENTRY(__epoll_pwait)


    movl    $__NR_epoll_pwait, %eax
    int     $0x80


    ret
END(__epoll_pwait)

