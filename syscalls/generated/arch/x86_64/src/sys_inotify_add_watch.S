#include <opa_asm.h>
    ENTRY(inotify_add_watch)


    movl    $__NR_inotify_add_watch, %eax
    syscall


    ret
END(inotify_add_watch)

