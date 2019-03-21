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
	.p2align 4,,15
	.globl	getSliceIndices
	.type	getSliceIndices, @function
getSliceIndices:
.LFB49:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movslq	%ecx, %rcx
	pushq	%rbx
	.cfi_def_cfa_offset 24
	.cfi_offset 3, -24
	leaq	(%rdx,%rcx,8), %rbx
	cmpq	$0, (%rbx)
	jle	.L22
	movslq	%r8d, %rax
	movslq	%r9d, %rbp
	xorl	%r10d, %r10d
	movslq	(%rsi,%rax,4), %r11
	.p2align 4,,10
	.p2align 3
.L16:
	testl	%r8d, %r8d
	movq	%r10, %rax
	jle	.L19
	xorl	%ecx, %ecx
	.p2align 4,,10
	.p2align 3
.L20:
	movslq	(%rsi,%rcx,4), %r9
	cqto
	addq	$1, %rcx
	idivq	%r9
	cmpl	%ecx, %r8d
	jg	.L20
.L19:
	cqto
	idivq	%r11
	cmpq	%rbp, %rdx
	je	.L24
.L18:
	addq	$1, %r10
	cmpq	%r10, (%rbx)
	jg	.L16
.L22:
	popq	%rbx
	.cfi_remember_state
	.cfi_def_cfa_offset 16
	popq	%rbp
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L24:
	.cfi_restore_state
	movq	%r10, (%rdi)
	addq	$8, %rdi
	jmp	.L18
	.cfi_endproc
.LFE49:
	.size	getSliceIndices, .-getSliceIndices
	.p2align 4,,15
	.globl	getFromInd
	.type	getFromInd, @function
getFromInd:
.LFB50:
	.cfi_startproc
	movslq	%ecx, %rcx
	xorl	%eax, %eax
	testq	%rcx, %rcx
	jle	.L25
	.p2align 4,,10
	.p2align 3
.L29:
	movq	(%rdx,%rax,8), %r8
	movsd	(%rsi,%r8,8), %xmm0
	movsd	%xmm0, (%rdi,%rax,8)
	addq	$1, %rax
	cmpq	%rcx, %rax
	jl	.L29
.L25:
	rep ret
	.cfi_endproc
.LFE50:
	.size	getFromInd, .-getFromInd
	.section	.rodata.str1.1
.LC9:
	.string	"Stack pointer (oracle): %p\n"
	.text
	.p2align 4,,15
	.globl	oracle2
	.type	oracle2, @function
oracle2:
.LFB51:
	.cfi_startproc
	subq	$24, %rsp
	.cfi_def_cfa_offset 32
	movl	$.LC9, %esi
	movl	$1, %edi
	leaq	15(%rsp), %rdx
	xorl	%eax, %eax
	call	__printf_chk
	leaq	15(%rsp), %rax
	addq	$24, %rsp
	.cfi_def_cfa_offset 8
	ret
	.cfi_endproc
.LFE51:
	.size	oracle2, .-oracle2
	.section	.rodata.str1.8,"aMS",@progbits,1
	.align 8
.LC10:
	.string	"GCC did not optimize this call.\n"
	.section	.rodata.str1.1
.LC11:
	.string	"YOLO"
	.text
	.p2align 4,,15
	.globl	getSliceInner
	.type	getSliceInner, @function
getSliceInner:
.LFB52:
	.cfi_startproc
	pushq	%r15
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
	xorl	%eax, %eax
	pushq	%r14
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	movq	%rdx, %r14
	pushq	%r13
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	movq	%rsi, %r13
	movl	$.LC9, %esi
	pushq	%r12
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	movq	%r9, %r12
	pushq	%rbp
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	movq	%r8, %rbp
	pushq	%rbx
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	movq	%rcx, %rbx
	subq	$40, %rsp
	.cfi_def_cfa_offset 96
	leaq	31(%rsp), %r15
	movq	%rdi, 8(%rsp)
	movl	$1, %edi
	movq	%r15, %rdx
	call	__printf_chk
	cmpq	%r12, %r15
	movq	8(%rsp), %r11
	jne	.L49
.L33:
	xorl	%eax, %eax
	movl	$.LC11, %esi
	movl	$1, %edi
	movq	%r11, 8(%rsp)
	call	__printf_chk
	movq	(%r14), %rcx
	movq	8(%rsp), %r11
	cmpq	%rbp, %rcx
	je	.L34
	movl	(%rbx), %eax
	testl	%eax, %eax
	jle	.L36
	leaq	-8(%r14), %rax
	leaq	-4(%rbx), %r15
	xorl	%r14d, %r14d
	movq	%rax, 8(%rsp)
	.p2align 4,,10
	.p2align 3
.L41:
	movq	8(%rsp), %rdx
	movq	%r11, %rdi
	movq	%r12, %r9
	movq	%rbp, %r8
	movq	%r15, %rcx
	movq	%r13, %rsi
	addl	$1, %r14d
	call	getSliceInner
	cmpl	%r14d, (%rbx)
	movq	%rax, %r11
	jg	.L41
.L36:
	addq	$40, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 56
	movq	%r11, %rax
	popq	%rbx
	.cfi_def_cfa_offset 48
	popq	%rbp
	.cfi_def_cfa_offset 40
	popq	%r12
	.cfi_def_cfa_offset 32
	popq	%r13
	.cfi_def_cfa_offset 24
	popq	%r14
	.cfi_def_cfa_offset 16
	popq	%r15
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L49:
	.cfi_restore_state
	testq	%r12, %r12
	je	.L33
	movl	$.LC10, %esi
	movl	$1, %edi
	xorl	%eax, %eax
	call	__printf_chk
	movq	8(%rsp), %r11
	jmp	.L33
	.p2align 4,,10
	.p2align 3
.L34:
	testq	%rbp, %rbp
	movq	0(%r13), %rdx
	jle	.L38
	salq	$3, %rbp
	xorl	%eax, %eax
	.p2align 4,,10
	.p2align 3
.L40:
	movsd	(%rdx,%rax), %xmm0
	movsd	%xmm0, (%r11,%rax)
	addq	$8, %rax
	cmpq	%rbp, %rax
	jne	.L40
	addq	%rax, %r11
	addq	%rax, %rdx
.L38:
	movl	(%rbx), %eax
	subl	$1, %eax
	cltq
	imulq	%rax, %rcx
	leaq	(%rdx,%rcx,8), %rax
	movq	%rax, 0(%r13)
	jmp	.L36
	.cfi_endproc
.LFE52:
	.size	getSliceInner, .-getSliceInner
	.p2align 4,,15
	.globl	getSlice
	.type	getSlice, @function
getSlice:
.LFB53:
	.cfi_startproc
	subq	$24, %rsp
	.cfi_def_cfa_offset 32
	movslq	%r9d, %r9
	movslq	%r8d, %r8
	movslq	32(%rsp), %rax
	movq	%rsi, 8(%rsp)
	leaq	-4(%rcx,%r8,4), %rcx
	movq	(%rdx,%r9,8), %rsi
	leaq	-8(%rdx,%r8,8), %rdx
	xorl	%r9d, %r9d
	imulq	%rsi, %rax
	movq	%rsi, %r8
	leaq	8(%rsp), %rsi
	salq	$3, %rax
	addq	%rax, 8(%rsp)
	call	getSliceInner
	addq	$24, %rsp
	.cfi_def_cfa_offset 8
	ret
	.cfi_endproc
.LFE53:
	.size	getSlice, .-getSlice
	.p2align 4,,15
	.globl	getSliceInner2
	.type	getSliceInner2, @function
getSliceInner2:
.LFB54:
	.cfi_startproc
	pushq	%r15
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
	movq	%rdi, %rax
	pushq	%r14
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	pushq	%r13
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	movq	%rsi, %r13
	pushq	%r12
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	movq	%r8, %r12
	pushq	%rbp
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	movq	%rcx, %rbp
	pushq	%rbx
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	subq	$8, %rsp
	.cfi_def_cfa_offset 64
	movq	(%rdx), %rdi
	cmpq	%r8, %rdi
	je	.L53
	movl	(%rcx), %ecx
	testl	%ecx, %ecx
	jle	.L55
	leaq	-4(%rbp), %r15
	leaq	-8(%rdx), %r14
	xorl	%ebx, %ebx
	.p2align 4,,10
	.p2align 3
.L60:
	movq	%r12, %r8
	movq	%r15, %rcx
	movq	%r14, %rdx
	movq	%r13, %rsi
	movq	%rax, %rdi
	addl	$1, %ebx
	call	getSliceInner2
	cmpl	%ebx, 0(%rbp)
	jg	.L60
.L55:
	addq	$8, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 56
	popq	%rbx
	.cfi_def_cfa_offset 48
	popq	%rbp
	.cfi_def_cfa_offset 40
	popq	%r12
	.cfi_def_cfa_offset 32
	popq	%r13
	.cfi_def_cfa_offset 24
	popq	%r14
	.cfi_def_cfa_offset 16
	popq	%r15
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L53:
	.cfi_restore_state
	testq	%rdi, %rdi
	movq	(%rsi), %rcx
	jle	.L57
	leaq	0(,%rdi,8), %rsi
	xorl	%edx, %edx
	.p2align 4,,10
	.p2align 3
.L59:
	movsd	(%rcx,%rdx), %xmm0
	movsd	%xmm0, (%rax,%rdx)
	addq	$8, %rdx
	cmpq	%rsi, %rdx
	jne	.L59
	addq	%rdx, %rax
	addq	%rdx, %rcx
.L57:
	movl	0(%rbp), %edx
	subl	$1, %edx
	movslq	%edx, %rdx
	imulq	%rdx, %rdi
	leaq	(%rcx,%rdi,8), %rdx
	movq	%rdx, 0(%r13)
	addq	$8, %rsp
	.cfi_def_cfa_offset 56
	popq	%rbx
	.cfi_def_cfa_offset 48
	popq	%rbp
	.cfi_def_cfa_offset 40
	popq	%r12
	.cfi_def_cfa_offset 32
	popq	%r13
	.cfi_def_cfa_offset 24
	popq	%r14
	.cfi_def_cfa_offset 16
	popq	%r15
	.cfi_def_cfa_offset 8
	ret
	.cfi_endproc
.LFE54:
	.size	getSliceInner2, .-getSliceInner2
	.p2align 4,,15
	.globl	getSlice2
	.type	getSlice2, @function
getSlice2:
.LFB55:
	.cfi_startproc
	subq	$24, %rsp
	.cfi_def_cfa_offset 32
	movslq	%r9d, %r9
	movslq	%r8d, %r8
	movslq	32(%rsp), %rax
	movq	%rsi, 8(%rsp)
	leaq	-4(%rcx,%r8,4), %rcx
	movq	(%rdx,%r9,8), %rsi
	leaq	-8(%rdx,%r8,8), %rdx
	imulq	%rsi, %rax
	movq	%rsi, %r8
	leaq	8(%rsp), %rsi
	salq	$3, %rax
	addq	%rax, 8(%rsp)
	call	getSliceInner2
	addq	$24, %rsp
	.cfi_def_cfa_offset 8
	ret
	.cfi_endproc
.LFE55:
	.size	getSlice2, .-getSlice2
	.p2align 4,,15
	.globl	factorial
	.type	factorial, @function
factorial:
.LFB56:
	.cfi_startproc
	pushq	%r12
	.cfi_def_cfa_offset 16
	.cfi_offset 12, -16
	xorl	%eax, %eax
	movq	%rsi, %r12
	movl	$.LC9, %esi
	pushq	%rbp
	.cfi_def_cfa_offset 24
	.cfi_offset 6, -24
	movl	%edi, %ebp
	movl	$1, %edi
	pushq	%rbx
	.cfi_def_cfa_offset 32
	.cfi_offset 3, -32
	subq	$16, %rsp
	.cfi_def_cfa_offset 48
	leaq	15(%rsp), %rbx
	movq	%rbx, %rdx
	call	__printf_chk
	testq	%r12, %r12
	je	.L65
	movslq	%ebx, %rax
	cmpq	%r12, %rax
	je	.L65
	movl	$.LC10, %esi
	movl	$1, %edi
	xorl	%eax, %eax
	call	__printf_chk
.L65:
	testl	%ebp, %ebp
	movl	$1, %eax
	je	.L66
	leal	-1(%rbp), %edi
	movslq	%ebx, %rsi
	call	factorial
	imull	%ebp, %eax
.L66:
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
.LFE56:
	.size	factorial, .-factorial
	.section	.rodata.str1.1
.LC12:
	.string	"Get indices once (naively)"
	.section	.rodata.str1.8
	.align 8
.LC13:
	.string	"Get slice using pre-stored indices 1000 times"
	.align 8
.LC14:
	.string	"Get slice directly using recursive scheme 1000 times"
	.section	.rodata.str1.1
.LC15:
	.string	"SUCCESS\n"
	.section	.rodata.str1.8
	.align 8
.LC16:
	.string	"ERROR: slice[%i]=%f!=%f=slice2[%i]\n"
	.section	.text.startup,"ax",@progbits
	.p2align 4,,15
	.globl	main
	.type	main, @function
main:
.LFB57:
	.cfi_startproc
	pushq	%r15
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
	movl	$16, %edi
	pushq	%r14
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	pushq	%r13
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	pushq	%r12
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	pushq	%rbp
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	pushq	%rbx
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	subq	$24, %rsp
	.cfi_def_cfa_offset 80
	call	malloc
	movl	$16, %edi
	movq	%rax, %rbp
	call	malloc
	movl	$40, %edi
	movq	%rax, %rbx
	call	malloc
	movl	$5, 0(%rbp)
	movq	%rax, %r13
	movl	$4, 4(%rbp)
	movl	$3, 8(%rbp)
	movl	$3, 12(%rbp)
	xorl	%edx, %edx
	movq	$1, (%rax)
.L75:
	movslq	0(%rbp,%rdx), %rcx
	imulq	0(%r13,%rdx,2), %rcx
	movq	%rcx, 8(%r13,%rdx,2)
	addq	$4, %rdx
	cmpq	$16, %rdx
	jne	.L75
	xorb	%dl, %dl
.L77:
	movl	0(%rbp,%rdx), %eax
	leal	-2(%rax), %ecx
	movl	%ecx, (%rbx,%rdx)
	addq	$4, %rdx
	cmpq	$16, %rdx
	jne	.L77
	movl	$32, %edi
	call	malloc
	xorl	%edx, %edx
	movq	%rax, %r11
.L82:
	movl	%edx, %esi
	movq	$1, (%r11,%rdx,8)
	xorl	%eax, %eax
.L80:
	cmpl	%eax, %esi
	je	.L78
	movslq	0(%rbp,%rax,4), %rcx
	imulq	(%r11,%rdx,8), %rcx
	movq	%rcx, (%r11,%rdx,8)
.L78:
	addq	$1, %rax
	cmpq	$4, %rax
	jne	.L80
	addq	$1, %rdx
	cmpq	$4, %rdx
	jne	.L82
	xorb	%al, %al
	xorl	%r14d, %r14d
.L84:
	movq	(%r11,%rax), %rdx
	cmpq	%rdx, %r14
	cmovl	%rdx, %r14
	addq	$8, %rax
	cmpq	$32, %rax
	jne	.L84
	movq	32(%r13), %r12
	movq	%r11, 8(%rsp)
	leaq	0(,%r12,8), %rdi
	call	malloc
	movq	%rax, %rbx
	xorl	%eax, %eax
	testq	%r12, %r12
	movq	8(%rsp), %r11
	jle	.L87
	.p2align 4,,10
	.p2align 3
.L104:
	cvtsi2sdq	%rax, %xmm0
	movsd	%xmm0, (%rbx,%rax,8)
	addq	$1, %rax
	cmpq	%r12, %rax
	jne	.L104
.L87:
	leaq	0(,%r14,8), %r15
	movq	%r11, 8(%rsp)
	movq	%r15, %rdi
	call	malloc
	movq	%r15, %rdi
	movq	%rax, %r14
	call	malloc
	xorl	%edi, %edi
	movq	%rax, %r12
	call	tprobe
	movq	%r15, %rdi
	call	malloc
	movl	$4, %ecx
	movq	%r13, %rdx
	movl	$2, %r9d
	movl	$1, %r8d
	movq	%rbp, %rsi
	movq	%rax, %rdi
	movq	%rax, %r15
	call	getSliceIndices
	movl	$.LC12, %edi
	call	tprobe
	xorl	%edi, %edi
	call	tprobe
	movq	8(%rsp), %r11
	movl	$1000, %ecx
	movq	8(%r11), %rax
	movq	%rax, 8(%rsp)
	movslq	%eax, %rdx
	.p2align 4,,10
	.p2align 3
.L86:
	xorl	%r8d, %r8d
	testq	%rdx, %rdx
	jle	.L92
	.p2align 4,,10
	.p2align 3
.L105:
	movq	(%r15,%r8,8), %rax
	movsd	(%rbx,%rax,8), %xmm0
	movsd	%xmm0, (%r14,%r8,8)
	addq	$1, %r8
	cmpq	%rdx, %r8
	jl	.L105
.L92:
	subl	$1, %ecx
	jne	.L86
	movl	$.LC13, %edi
	movl	$1000, %r15d
	call	tprobe
	xorl	%edi, %edi
	call	tprobe
	.p2align 4,,10
	.p2align 3
.L95:
	movl	$2, (%rsp)
	movl	$1, %r9d
	movl	$4, %r8d
	movq	%rbp, %rcx
	movq	%r13, %rdx
	movq	%rbx, %rsi
	movq	%r12, %rdi
	call	getSlice
	subl	$1, %r15d
	jne	.L95
	movl	$.LC14, %edi
	xorl	%ebx, %ebx
	call	tprobe
	xorl	%esi, %esi
	movl	$20, %edi
	call	factorial
	cmpq	$0, 8(%rsp)
	movl	$1, %eax
	jg	.L106
	jmp	.L102
	.p2align 4,,10
	.p2align 3
.L110:
	jne	.L103
	addq	$1, %rbx
	cmpq	8(%rsp), %rbx
	.p2align 4,,2
	je	.L109
.L106:
	movsd	(%r14,%rbx,8), %xmm0
	movl	%ebx, %edx
	movsd	(%r12,%rbx,8), %xmm1
	ucomisd	%xmm1, %xmm0
	jnp	.L110
.L103:
	movl	%edx, %ecx
	movl	$.LC16, %esi
	movl	$1, %edi
	movl	$2, %eax
	addq	$1, %rbx
	call	__printf_chk
	xorl	%eax, %eax
	cmpq	8(%rsp), %rbx
	jne	.L106
.L109:
	testl	%eax, %eax
	jne	.L102
.L97:
	addq	$24, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 56
	xorl	%eax, %eax
	popq	%rbx
	.cfi_def_cfa_offset 48
	popq	%rbp
	.cfi_def_cfa_offset 40
	popq	%r12
	.cfi_def_cfa_offset 32
	popq	%r13
	.cfi_def_cfa_offset 24
	popq	%r14
	.cfi_def_cfa_offset 16
	popq	%r15
	.cfi_def_cfa_offset 8
	ret
.L102:
	.cfi_restore_state
	movl	$.LC15, %esi
	movl	$1, %edi
	xorl	%eax, %eax
	call	__printf_chk
	jmp	.L97
	.cfi_endproc
.LFE57:
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
