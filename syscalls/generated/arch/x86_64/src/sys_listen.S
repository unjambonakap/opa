#include <opa_asm.h>
    ENTRY(listen)


    movl    $__NR_listen, %eax
    syscall


    ret
END(listen)

