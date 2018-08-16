#include <opa_asm.h>
    ENTRY(__fcntl64)


    movl    $__NR_fcntl64, %eax
    int     $0x80


    ret
END(__fcntl64)

