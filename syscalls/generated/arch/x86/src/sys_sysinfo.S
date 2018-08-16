#include <opa_asm.h>
    ENTRY(sysinfo)


    movl    $__NR_sysinfo, %eax
    int     $0x80


    ret
END(sysinfo)

