#include <opa_asm.h>
    ENTRY(inotify_rm_watch)


    movl    $__NR_inotify_rm_watch, %eax
    int     $0x80


    ret
END(inotify_rm_watch)

