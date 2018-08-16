#include <opa_asm.h>
    ENTRY(setrlimit)


    movl    $__NR_setrlimit, %eax
    syscall


    ret
END(setrlimit)

