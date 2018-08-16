#include <opa_asm.h>
    ENTRY(__set_tid_address)


    movl    $__NR_set_tid_address, %eax
    int     $0x80


    ret
END(__set_tid_address)

