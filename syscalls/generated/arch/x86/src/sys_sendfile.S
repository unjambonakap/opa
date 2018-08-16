#include <opa_asm.h>
    ENTRY(sendfile)


    movl    $__NR_sendfile, %eax
    int     $0x80


    ret
END(sendfile)

