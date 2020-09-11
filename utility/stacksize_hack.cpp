
void main2(){
    // code here
}

static void run_with_stack_size(void (*func)(void), size_t stsize) {
    char *stack, *send;
    stack=(char *)malloc(stsize);
    send=stack+stsize-16;
    send=(char *)((uintptr_t)send/16*16);
    asm volatile(
        "mov %%rsp, (%0)\n"
        "mov %0, %%rsp\n"
        :
        : "r" (send));
    func();
    asm volatile(
         "mov (%0), %%rsp\n"
         :
         : "r" (send));
    free(stack);
}
int main() {
    run_with_stack_size(main2, 1024*1024*1024);
    return 0;
}
