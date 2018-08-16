#define OPA_SYMNAME(name) opa_KernSym_ ## name 
OPA_SYMVAL(start_thread)
OPA_SYMVAL(start_thread_ia32)
OPA_SYMVAL(kernel_thread)
OPA_SYMVAL(access_remote_vm)
OPA_SYMVAL(do_mmap_pgoff)
OPA_SYMVAL(__mm_populate)
#undef OPA_SYNAME
