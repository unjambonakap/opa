#include <opa_asm.h>
    ENTRY(sysinfo)


    movl    $__NR_sysinfo, %eax
    syscall


    ret
END(sysinfo)

