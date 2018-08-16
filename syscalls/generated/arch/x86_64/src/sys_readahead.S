#include <opa_asm.h>
    ENTRY(readahead)


    movl    $__NR_readahead, %eax
    syscall


    ret
END(readahead)

