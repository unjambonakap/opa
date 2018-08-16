#include <opa_asm.h>
    ENTRY(fdatasync)


    movl    $__NR_fdatasync, %eax
    int     $0x80


    ret
END(fdatasync)

