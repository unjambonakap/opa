#include <opa_asm.h>
    ENTRY(epoll_ctl)

mov %rcx, %r10

    movl    $__NR_epoll_ctl, %eax
    syscall


    ret
END(epoll_ctl)

