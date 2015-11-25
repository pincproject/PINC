	.file	"main.c"
	.section	.rodata.str1.1,"aMS",@progbits,1
.LC1:
	.string	"Time: "
.LC2:
	.string	"%.2fs"
.LC4:
	.string	"%.1fms"
.LC6:
	.string	"%.1fus"
.LC7:
	.string	"%ldns"
.LC8:
	.string	" %s\n"
	.text
	.p2align 4,,15
	.globl	tprobe
	.type	tprobe, @function
tprobe:
.LFB48:
	.cfi_startproc
	pushq	%r12
	.cfi_def_cfa_offset 16
	.cfi_offset 12, -16
	movq	%rdi, %r12
	movl	$2, %edi
	pushq	%rbp
	.cfi_def_cfa_offset 24
	.cfi_offset 6, -24
	pushq	%rbx
	.cfi_def_cfa_offset 32
	.cfi_offset 3, -32
	subq	$16, %rsp
	.cfi_def_cfa_offset 48
	movq	%rsp, %rsi
	call	clock_gettime
	testq	%r12, %r12
	je	.L11
	movq	(%rsp), %rbx
	movq	8(%rsp), %rbp
	subq	previous.3046(%rip), %rbx
	subq	previous.3046+8(%rip), %rbp
	js	.L12
.L4:
	xorl	%eax, %eax
	movl	$.LC1, %esi
	movl	$1, %edi
	call	__printf_chk
	testq	%rbx, %rbx
	jle	.L5
	cvtsi2sdq	%rbp, %xmm1
	movl	$.LC2, %esi
	cvtsi2sdq	%rbx, %xmm0
	movl	$1, %edi
	movl	$1, %eax
	divsd	.LC0(%rip), %xmm1
	addsd	%xmm1, %xmm0
	call	__printf_chk
.L6:
	movq	%r12, %rdx
	movl	$.LC8, %esi
	movl	$1, %edi
	xorl	%eax, %eax
	call	__printf_chk
	addq	$16, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 32
	popq	%rbx
	.cfi_def_cfa_offset 24
	popq	%rbp
	.cfi_def_cfa_offset 16
	popq	%r12
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L5:
	.cfi_restore_state
	cmpq	$1000000, %rbp
	jg	.L13
	cmpq	$1000, %rbp
	jle	.L8
	cvtsi2sdq	%rbp, %xmm0
	movl	$.LC6, %esi
	movl	$1, %edi
	movl	$1, %eax
	divsd	.LC5(%rip), %xmm0
	call	__printf_chk
	jmp	.L6
	.p2align 4,,10
	.p2align 3
.L12:
	cvtsi2sdq	%rbp, %xmm0
	subq	$1, %rbx
	addsd	.LC0(%rip), %xmm0
	cvttsd2siq	%xmm0, %rbp
	jmp	.L4
	.p2align 4,,10
	.p2align 3
.L8:
	movq	%rbp, %rdx
	movl	$.LC7, %esi
	movl	$1, %edi
	xorl	%eax, %eax
	call	__printf_chk
	jmp	.L6
	.p2align 4,,10
	.p2align 3
.L13:
	cvtsi2sdq	%rbp, %xmm0
	movl	$.LC4, %esi
	movl	$1, %edi
	movl	$1, %eax
	divsd	.LC3(%rip), %xmm0
	call	__printf_chk
	jmp	.L6
	.p2align 4,,10
	.p2align 3
.L11:
	movq	(%rsp), %rax
	movq	8(%rsp), %rdx
	movq	%rax, previous.3046(%rip)
	movq	%rdx, previous.3046+8(%rip)
	addq	$16, %rsp
	.cfi_def_cfa_offset 32
	popq	%rbx
	.cfi_def_cfa_offset 24
	popq	%rbp
	.cfi_def_cfa_offset 16
	popq	%r12
	.cfi_def_cfa_offset 8
	ret
	.cfi_endproc
.LFE48:
	.size	tprobe, .-tprobe
	.section	.rodata.str1.1
.LC9:
	.string	"Stack pointer: %p\n"
	.text
	.p2align 4,,15
	.globl	printSP2
	.type	printSP2, @function
printSP2:
.LFB49:
	.cfi_startproc
	subq	$24, %rsp
	.cfi_def_cfa_offset 32
	movl	$.LC9, %esi
	movl	$1, %edi
	movq	%rsp, %rdx
	xorl	%eax, %eax
	movq	$0, (%rsp)
	call	__printf_chk
	addq	$24, %rsp
	.cfi_def_cfa_offset 8
	ret
	.cfi_endproc
.LFE49:
	.size	printSP2, .-printSP2
	.p2align 4,,15
	.globl	printSP
	.type	printSP, @function
printSP:
.LFB50:
	.cfi_startproc
	movq	%rsp, %rdx
	movl	$.LC9, %esi
	movl	$1, %edi
	xorl	%eax, %eax
	jmp	__printf_chk
	.cfi_endproc
.LFE50:
	.size	printSP, .-printSP
	.section	.rodata.str1.1
.LC10:
	.string	"Frame address: %p\n"
	.text
	.p2align 4,,15
	.globl	factorial
	.type	factorial, @function
factorial:
.LFB51:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%r13
	pushq	%r12
	pushq	%rbx
	.cfi_offset 13, -24
	.cfi_offset 12, -32
	.cfi_offset 3, -40
	movq	%rdi, %rbx
	subq	$24, %rsp
	cmpq	$5, %rdi
	je	.L33
	xorl	%eax, %eax
	movq	%rbp, %rdx
	movl	$.LC10, %esi
	movl	$1, %edi
	call	__printf_chk
	testq	%rbx, %rbx
	movl	$1, %eax
	jne	.L34
	addq	$24, %rsp
	popq	%rbx
	popq	%r12
	popq	%r13
	popq	%rbp
	.cfi_remember_state
	.cfi_def_cfa 7, 8
	ret
	.p2align 4,,10
	.p2align 3
.L33:
	.cfi_restore_state
	movq	%rbp, %rdx
	movl	$.LC10, %esi
	movl	$1, %edi
	xorl	%eax, %eax
	movl	$8, -44(%rbp)
	movl	$3, %r13d
	call	__printf_chk
	movl	$4, %r12d
	movq	%rbp, %rdx
	movl	$.LC10, %esi
	movl	$1, %edi
	xorl	%eax, %eax
	call	__printf_chk
.L19:
	xorl	%eax, %eax
	movq	%rbp, %rdx
	movl	$.LC10, %esi
	movl	$1, %edi
	call	__printf_chk
	testq	%r13, %r13
	movl	$1, %eax
	jne	.L24
.L25:
	imulq	%r12, %rax
.L23:
	addq	$24, %rsp
	imulq	%rbx, %rax
	popq	%rbx
	popq	%r12
	popq	%r13
	popq	%rbp
	.cfi_remember_state
	.cfi_def_cfa 7, 8
	ret
	.p2align 4,,10
	.p2align 3
.L34:
	.cfi_restore_state
	leaq	-1(%rbx), %r12
	cmpq	$5, %r12
	jne	.L21
	movq	%rbp, %rdx
	movl	$.LC10, %esi
	movl	$1, %edi
	xorl	%eax, %eax
	movl	$8, -40(%rbp)
	call	__printf_chk
.L22:
	leaq	-2(%rbx), %r13
	cmpq	$5, %r13
	jne	.L19
	movq	%rbp, %rdx
	movl	$.LC10, %esi
	movl	$1, %edi
	xorl	%eax, %eax
	movl	$8, -36(%rbp)
	call	__printf_chk
.L24:
	leaq	-3(%rbx), %rdi
	call	factorial
	imulq	%r13, %rax
	jmp	.L25
.L21:
	xorl	%eax, %eax
	movq	%rbp, %rdx
	movl	$.LC10, %esi
	movl	$1, %edi
	call	__printf_chk
	testq	%r12, %r12
	movl	$1, %eax
	je	.L23
	jmp	.L22
	.cfi_endproc
.LFE51:
	.size	factorial, .-factorial
	.section	.rodata.str1.1
.LC11:
	.string	"%i!=%li\n"
	.section	.text.startup,"ax",@progbits
	.p2align 4,,15
	.globl	main
	.type	main, @function
main:
.LFB52:
	.cfi_startproc
	subq	$8, %rsp
	.cfi_def_cfa_offset 16
	movl	$10, %edi
	call	factorial
	movl	$10, %edx
	movq	%rax, %rcx
	movl	$.LC11, %esi
	movl	$1, %edi
	xorl	%eax, %eax
	call	__printf_chk
	xorl	%eax, %eax
	addq	$8, %rsp
	.cfi_def_cfa_offset 8
	ret
	.cfi_endproc
.LFE52:
	.size	main, .-main
	.local	previous.3046
	.comm	previous.3046,16,16
	.section	.rodata.cst8,"aM",@progbits,8
	.align 8
.LC0:
	.long	0
	.long	1104006501
	.align 8
.LC3:
	.long	0
	.long	1093567616
	.align 8
.LC5:
	.long	0
	.long	1083129856
	.ident	"GCC: (Ubuntu 4.8.2-19ubuntu1) 4.8.2"
	.section	.note.GNU-stack,"",@progbits
