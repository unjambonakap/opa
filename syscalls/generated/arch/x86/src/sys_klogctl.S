#include <opa_asm.h>
    ENTRY(klogctl)


    movl    $__NR_syslog, %eax
    int     $0x80


    ret
END(klogctl)

