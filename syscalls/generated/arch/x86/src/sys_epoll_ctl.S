#include <opa_asm.h>
    ENTRY(epoll_ctl)


    movl    $__NR_epoll_ctl, %eax
    int     $0x80


    ret
END(epoll_ctl)

