#include <opa_asm.h>
    ENTRY(fstatat64)


    movl    $__NR_fstatat64, %eax
    int     $0x80


    ret
END(fstatat64)

