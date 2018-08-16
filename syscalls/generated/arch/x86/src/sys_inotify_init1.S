#include <opa_asm.h>
    ENTRY(inotify_init1)


    movl    $__NR_inotify_init1, %eax
    int     $0x80


    ret
END(inotify_init1)

