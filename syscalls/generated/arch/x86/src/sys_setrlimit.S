#include <opa_asm.h>
    ENTRY(setrlimit)


    movl    $__NR_setrlimit, %eax
    int     $0x80


    ret
END(setrlimit)

