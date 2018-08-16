#include <opa_asm.h>
    ENTRY(sethostname)


    movl    $__NR_sethostname, %eax
    int     $0x80


    ret
END(sethostname)

