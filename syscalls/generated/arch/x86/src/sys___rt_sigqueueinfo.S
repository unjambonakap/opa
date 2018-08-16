#include <opa_asm.h>
    ENTRY(__rt_sigqueueinfo)


    movl    $__NR_rt_sigqueueinfo, %eax
    int     $0x80


    ret
END(__rt_sigqueueinfo)

