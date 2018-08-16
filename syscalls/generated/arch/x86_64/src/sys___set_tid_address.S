#include <opa_asm.h>
    ENTRY(__set_tid_address)


    movl    $__NR_set_tid_address, %eax
    syscall


    ret
END(__set_tid_address)

