#include <opa_asm.h>
    ENTRY(sendfile)

mov %rcx, %r10

    movl    $__NR_sendfile, %eax
    syscall


    ret
END(sendfile)

