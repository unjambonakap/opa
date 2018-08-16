#include <opa_asm.h>
    ENTRY(fsync)


    movl    $__NR_fsync, %eax
    int     $0x80


    ret
END(fsync)

