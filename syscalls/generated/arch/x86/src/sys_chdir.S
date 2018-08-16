#include <opa_asm.h>
    ENTRY(chdir)


    movl    $__NR_chdir, %eax
    int     $0x80


    ret
END(chdir)

