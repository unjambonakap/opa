#include <opa_asm.h>
    ENTRY(sendfile64)


    movl    $__NR_sendfile64, %eax
    int     $0x80


    ret
END(sendfile64)

