#include <opa_asm.h>
    ENTRY(__set_thread_area)


    movl    $__NR_set_thread_area, %eax
    int     $0x80


    ret
END(__set_thread_area)

