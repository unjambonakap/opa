#include <opa_asm.h>
    ENTRY(munmap)


    movl    $__NR_munmap, %eax
    syscall


    ret
END(munmap)

