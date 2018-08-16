#include <opa_asm.h>
    ENTRY(_exit)


    movl    $__NR_exit_group, %eax
    syscall


    ret
END(_exit)

