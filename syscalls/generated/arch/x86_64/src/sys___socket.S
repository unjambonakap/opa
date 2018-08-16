#include <opa_asm.h>
    ENTRY(__socket)


    movl    $__NR_socket, %eax
    syscall


    ret
END(__socket)

