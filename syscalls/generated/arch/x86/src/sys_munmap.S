#include <opa_asm.h>
    ENTRY(munmap)


    movl    $__NR_munmap, %eax
    int     $0x80


    ret
END(munmap)

