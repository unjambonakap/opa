#include <opa_asm.h>
    ENTRY(klogctl)


    movl    $__NR_syslog, %eax
    syscall


    ret
END(klogctl)

