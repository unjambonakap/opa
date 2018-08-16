#include <opa_asm.h>
    ENTRY(fchdir)


    movl    $__NR_fchdir, %eax
    int     $0x80


    ret
END(fchdir)

