#include <opa_asm.h>
    ENTRY(inotify_init1)


    movl    $__NR_inotify_init1, %eax
    syscall


    ret
END(inotify_init1)

