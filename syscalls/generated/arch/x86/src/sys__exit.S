#include <opa_asm.h>
    ENTRY(_exit)


    movl    $__NR_exit_group, %eax
    int     $0x80


    ret
END(_exit)

