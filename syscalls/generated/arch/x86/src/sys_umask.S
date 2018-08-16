#include <opa_asm.h>
    ENTRY(umask)


    movl    $__NR_umask, %eax
    int     $0x80


    ret
END(umask)

