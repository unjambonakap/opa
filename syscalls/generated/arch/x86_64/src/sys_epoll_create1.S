#include <opa_asm.h>
    ENTRY(epoll_create1)


    movl    $__NR_epoll_create1, %eax
    syscall


    ret
END(epoll_create1)

