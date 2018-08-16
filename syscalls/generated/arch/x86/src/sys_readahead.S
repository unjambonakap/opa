#include <opa_asm.h>
    ENTRY(readahead)


    movl    $__NR_readahead, %eax
    int     $0x80


    ret
END(readahead)

