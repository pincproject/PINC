	.file	"halo.c"
	.text
	.p2align 4,,15
	.globl	getSlice
	.type	getSlice, @function
getSlice:
.LFB0:
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
	pushq	%r12
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	pushq	%rbp
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	movq	%rcx, %rbp
	pushq	%rbx
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	movq	%rsi, %rbx
	subq	$56, %rsp
	.cfi_def_cfa_offset 112
	movq	(%rdx), %rcx
	movq	%rdx, 8(%rsp)
	testq	%rcx, %rcx
	jle	.L2
	movq	(%rsi), %rsi
	leaq	0(,%rcx,8), %rdi
	xorl	%edx, %edx
.L4:
	movsd	(%rsi,%rdx), %xmm0
	movsd	%xmm0, (%rax,%rdx)
	addq	$8, %rdx
	cmpq	%rdi, %rdx
	jne	.L4
	addq	%rdx, %rax
	addq	%rsi, %rdx
	cmpq	$1, %rcx
	movq	%rdx, (%rbx)
	je	.L5
.L2:
	movl	0(%rbp), %r12d
	testl	%r12d, %r12d
	jle	.L101
	movq	8(%rsp), %rdi
	leaq	-36(%rbp), %r13
	movl	$0, 24(%rsp)
	movq	-8(%rdi), %rcx
	leaq	-72(%rdi), %r12
.L7:
	testq	%rcx, %rcx
	jle	.L95
	movq	(%rbx), %rsi
	leaq	0(,%rcx,8), %rdi
	xorl	%edx, %edx
.L96:
	movsd	(%rsi,%rdx), %xmm0
	movsd	%xmm0, (%rax,%rdx)
	addq	$8, %rdx
	cmpq	%rdi, %rdx
	jne	.L96
	addq	%rdx, %rax
	addq	%rsi, %rdx
	cmpq	$1, %rcx
	movq	%rdx, (%rbx)
	je	.L9
.L95:
	movl	-4(%rbp), %r11d
	testl	%r11d, %r11d
	jle	.L10
	movq	8(%rsp), %rdi
	movl	$0, 28(%rsp)
	movq	-16(%rdi), %rcx
.L11:
	testq	%rcx, %rcx
	jle	.L88
	movq	(%rbx), %rsi
	leaq	0(,%rcx,8), %rdi
	xorl	%edx, %edx
.L89:
	movsd	(%rsi,%rdx), %xmm0
	movsd	%xmm0, (%rax,%rdx)
	addq	$8, %rdx
	cmpq	%rdi, %rdx
	jne	.L89
	addq	%rdx, %rax
	addq	%rsi, %rdx
	cmpq	$1, %rcx
	movq	%rdx, (%rbx)
	je	.L13
.L88:
	movl	-8(%rbp), %r10d
	testl	%r10d, %r10d
	jle	.L14
	movq	8(%rsp), %rdi
	movl	$0, 32(%rsp)
	movq	-24(%rdi), %rcx
.L15:
	testq	%rcx, %rcx
	jle	.L81
	movq	(%rbx), %rsi
	leaq	0(,%rcx,8), %rdi
	xorl	%edx, %edx
.L82:
	movsd	(%rsi,%rdx), %xmm0
	movsd	%xmm0, (%rax,%rdx)
	addq	$8, %rdx
	cmpq	%rdi, %rdx
	jne	.L82
	addq	%rdx, %rax
	addq	%rsi, %rdx
	cmpq	$1, %rcx
	movq	%rdx, (%rbx)
	je	.L17
.L81:
	movl	-12(%rbp), %r9d
	testl	%r9d, %r9d
	jle	.L18
	movq	8(%rsp), %rdi
	movl	$0, 36(%rsp)
	movq	-32(%rdi), %rcx
.L19:
	testq	%rcx, %rcx
	jle	.L74
	movq	(%rbx), %rsi
	leaq	0(,%rcx,8), %rdi
	xorl	%edx, %edx
.L75:
	movsd	(%rsi,%rdx), %xmm0
	movsd	%xmm0, (%rax,%rdx)
	addq	$8, %rdx
	cmpq	%rdi, %rdx
	jne	.L75
	addq	%rdx, %rax
	addq	%rsi, %rdx
	cmpq	$1, %rcx
	movq	%rdx, (%rbx)
	je	.L21
.L74:
	movl	-16(%rbp), %r8d
	testl	%r8d, %r8d
	jle	.L22
	movq	8(%rsp), %rdi
	movl	$0, 40(%rsp)
	movq	-40(%rdi), %rcx
.L23:
	testq	%rcx, %rcx
	jle	.L67
	movq	(%rbx), %rsi
	leaq	0(,%rcx,8), %rdi
	xorl	%edx, %edx
.L68:
	movsd	(%rsi,%rdx), %xmm0
	movsd	%xmm0, (%rax,%rdx)
	addq	$8, %rdx
	cmpq	%rdi, %rdx
	jne	.L68
	addq	%rdx, %rax
	addq	%rsi, %rdx
	cmpq	$1, %rcx
	movq	%rdx, (%rbx)
	je	.L25
.L67:
	movl	-20(%rbp), %edi
	testl	%edi, %edi
	jle	.L26
	movq	8(%rsp), %rdi
	movl	$0, 44(%rsp)
	movq	-48(%rdi), %rdi
.L27:
	testq	%rdi, %rdi
	jle	.L60
	movq	(%rbx), %rcx
	leaq	0(,%rdi,8), %rsi
	xorl	%edx, %edx
	.p2align 4,,10
	.p2align 3
.L61:
	movsd	(%rcx,%rdx), %xmm0
	movsd	%xmm0, (%rax,%rdx)
	addq	$8, %rdx
	cmpq	%rsi, %rdx
	jne	.L61
	addq	%rdx, %rax
	addq	%rcx, %rdx
	cmpq	$1, %rdi
	movq	%rdx, (%rbx)
	je	.L29
.L60:
	movl	-24(%rbp), %esi
	testl	%esi, %esi
	jle	.L30
	movq	8(%rsp), %rdi
	movl	$0, 20(%rsp)
	movq	-56(%rdi), %rdi
	.p2align 4,,10
	.p2align 3
.L31:
	testq	%rdi, %rdi
	jle	.L53
	movq	(%rbx), %rcx
	leaq	0(,%rdi,8), %rsi
	xorl	%edx, %edx
	.p2align 4,,10
	.p2align 3
.L54:
	movsd	(%rcx,%rdx), %xmm0
	movsd	%xmm0, (%rax,%rdx)
	addq	$8, %rdx
	cmpq	%rsi, %rdx
	jne	.L54
	addq	%rdx, %rax
	addq	%rcx, %rdx
	cmpq	$1, %rdi
	movq	%rdx, (%rbx)
	je	.L33
.L53:
	movl	-28(%rbp), %ecx
	testl	%ecx, %ecx
	jle	.L34
	movq	8(%rsp), %rdi
	xorl	%r14d, %r14d
	movq	-64(%rdi), %rdi
	.p2align 4,,10
	.p2align 3
.L35:
	testq	%rdi, %rdi
	jle	.L46
	movq	(%rbx), %rcx
	leaq	0(,%rdi,8), %rsi
	xorl	%edx, %edx
	.p2align 4,,10
	.p2align 3
.L47:
	movsd	(%rcx,%rdx), %xmm0
	movsd	%xmm0, (%rax,%rdx)
	addq	$8, %rdx
	cmpq	%rsi, %rdx
	jne	.L47
	addq	%rdx, %rax
	addq	%rcx, %rdx
	cmpq	$1, %rdi
	movq	%rdx, (%rbx)
	je	.L37
.L46:
	movl	-32(%rbp), %edx
	testl	%edx, %edx
	jle	.L38
	xorl	%r15d, %r15d
	.p2align 4,,10
	.p2align 3
.L40:
	movq	%r13, %rcx
	movq	%r12, %rdx
	movq	%rbx, %rsi
	movq	%rax, %rdi
	addl	$1, %r15d
	call	getSlice
	cmpl	-32(%rbp), %r15d
	jl	.L40
	movq	8(%rsp), %rdi
	movq	-64(%rdi), %rdi
.L38:
	testq	%rdi, %rdi
	jle	.L42
.L41:
	movq	(%rbx), %rcx
	leaq	0(,%rdi,8), %rsi
	xorl	%edx, %edx
	.p2align 4,,10
	.p2align 3
.L44:
	movsd	(%rcx,%rdx), %xmm0
	movsd	%xmm0, (%rax,%rdx)
	addq	$8, %rdx
	cmpq	%rsi, %rdx
	jne	.L44
	addq	%rdx, %rax
	addq	%rcx, %rdx
	movq	%rdx, (%rbx)
.L42:
	addl	$1, %r14d
	cmpl	-28(%rbp), %r14d
	jl	.L35
	movq	8(%rsp), %rdi
	movq	-56(%rdi), %rdi
.L34:
	testq	%rdi, %rdi
	jle	.L49
.L48:
	movq	(%rbx), %rcx
	leaq	0(,%rdi,8), %rsi
	xorl	%edx, %edx
	.p2align 4,,10
	.p2align 3
.L51:
	movsd	(%rcx,%rdx), %xmm0
	movsd	%xmm0, (%rax,%rdx)
	addq	$8, %rdx
	cmpq	%rsi, %rdx
	jne	.L51
	addq	%rdx, %rax
	addq	%rcx, %rdx
	movq	%rdx, (%rbx)
.L49:
	addl	$1, 20(%rsp)
	movl	20(%rsp), %esi
	cmpl	-24(%rbp), %esi
	jl	.L31
	movq	8(%rsp), %rdi
	movq	-48(%rdi), %rdi
.L30:
	testq	%rdi, %rdi
	jle	.L56
.L55:
	movq	(%rbx), %rcx
	leaq	0(,%rdi,8), %rsi
	xorl	%edx, %edx
	.p2align 4,,10
	.p2align 3
.L58:
	movsd	(%rcx,%rdx), %xmm0
	movsd	%xmm0, (%rax,%rdx)
	addq	$8, %rdx
	cmpq	%rsi, %rdx
	jne	.L58
	addq	%rdx, %rax
	addq	%rcx, %rdx
	movq	%rdx, (%rbx)
.L56:
	addl	$1, 44(%rsp)
	movl	44(%rsp), %esi
	cmpl	-20(%rbp), %esi
	jl	.L27
	movq	8(%rsp), %rdi
	movq	-40(%rdi), %rcx
.L26:
	testq	%rcx, %rcx
	jle	.L63
.L62:
	movq	(%rbx), %rsi
	leaq	0(,%rcx,8), %rdi
	xorl	%edx, %edx
.L65:
	movsd	(%rsi,%rdx), %xmm0
	movsd	%xmm0, (%rax,%rdx)
	addq	$8, %rdx
	cmpq	%rdi, %rdx
	jne	.L65
	addq	%rdx, %rax
	addq	%rsi, %rdx
	movq	%rdx, (%rbx)
.L63:
	addl	$1, 40(%rsp)
	movl	40(%rsp), %edi
	cmpl	-16(%rbp), %edi
	jl	.L23
	movq	8(%rsp), %rdi
	movq	-32(%rdi), %rcx
.L22:
	testq	%rcx, %rcx
	jle	.L70
.L69:
	movq	(%rbx), %rsi
	leaq	0(,%rcx,8), %rdi
	xorl	%edx, %edx
.L72:
	movsd	(%rsi,%rdx), %xmm0
	movsd	%xmm0, (%rax,%rdx)
	addq	$8, %rdx
	cmpq	%rdi, %rdx
	jne	.L72
	addq	%rdx, %rax
	addq	%rsi, %rdx
	movq	%rdx, (%rbx)
.L70:
	addl	$1, 36(%rsp)
	movl	36(%rsp), %edi
	cmpl	-12(%rbp), %edi
	jl	.L19
	movq	8(%rsp), %rdi
	movq	-24(%rdi), %rcx
.L18:
	testq	%rcx, %rcx
	jle	.L77
.L76:
	movq	(%rbx), %rsi
	leaq	0(,%rcx,8), %rdi
	xorl	%edx, %edx
.L79:
	movsd	(%rsi,%rdx), %xmm0
	movsd	%xmm0, (%rax,%rdx)
	addq	$8, %rdx
	cmpq	%rdi, %rdx
	jne	.L79
	addq	%rdx, %rax
	addq	%rsi, %rdx
	movq	%rdx, (%rbx)
.L77:
	addl	$1, 32(%rsp)
	movl	32(%rsp), %edi
	cmpl	-8(%rbp), %edi
	jl	.L15
	movq	8(%rsp), %rdi
	movq	-16(%rdi), %rcx
.L14:
	testq	%rcx, %rcx
	jle	.L84
.L83:
	movq	(%rbx), %rsi
	leaq	0(,%rcx,8), %rdi
	xorl	%edx, %edx
.L86:
	movsd	(%rsi,%rdx), %xmm0
	movsd	%xmm0, (%rax,%rdx)
	addq	$8, %rdx
	cmpq	%rdi, %rdx
	jne	.L86
	addq	%rdx, %rax
	addq	%rsi, %rdx
	movq	%rdx, (%rbx)
.L84:
	addl	$1, 28(%rsp)
	movl	28(%rsp), %edi
	cmpl	-4(%rbp), %edi
	jl	.L11
	movq	8(%rsp), %rdi
	movq	-8(%rdi), %rcx
.L10:
	testq	%rcx, %rcx
	jg	.L90
.L91:
	addl	$1, 24(%rsp)
	movl	24(%rsp), %edi
	cmpl	%edi, 0(%rbp)
	jg	.L7
	movq	8(%rsp), %rdi
	movq	(%rdi), %rdx
	jmp	.L6
.L5:
	movslq	0(%rbp), %rdx
	salq	$3, %rdx
	addq	%rdx, (%rbx)
	movl	$1, %edx
.L97:
	movq	(%rbx), %rcx
	leaq	0(,%rdx,8), %rsi
	xorl	%edx, %edx
.L100:
	movsd	(%rcx,%rdx), %xmm0
	movsd	%xmm0, (%rax,%rdx)
	addq	$8, %rdx
	cmpq	%rsi, %rdx
	jne	.L100
	addq	%rdx, %rax
	addq	%rcx, %rdx
	movq	%rdx, (%rbx)
.L98:
	addq	$56, %rsp
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
.L37:
	.cfi_restore_state
	movslq	-32(%rbp), %rdx
	salq	$3, %rdx
	addq	%rdx, (%rbx)
	jmp	.L41
.L33:
	movslq	-28(%rbp), %rdx
	salq	$3, %rdx
	addq	%rdx, (%rbx)
	jmp	.L48
.L29:
	movslq	-24(%rbp), %rdx
	salq	$3, %rdx
	addq	%rdx, (%rbx)
	jmp	.L55
.L25:
	movslq	-20(%rbp), %rdx
	salq	$3, %rdx
	addq	%rdx, (%rbx)
	jmp	.L62
.L101:
	movq	%rcx, %rdx
.L6:
	testq	%rdx, %rdx
	jg	.L97
	jmp	.L98
.L9:
	movslq	-4(%rbp), %rdx
	salq	$3, %rdx
	addq	%rdx, (%rbx)
.L90:
	movq	(%rbx), %rsi
	leaq	0(,%rcx,8), %rdi
	xorl	%edx, %edx
.L93:
	movsd	(%rsi,%rdx), %xmm0
	movsd	%xmm0, (%rax,%rdx)
	addq	$8, %rdx
	cmpq	%rdi, %rdx
	jne	.L93
	addq	%rdx, %rax
	addq	%rsi, %rdx
	movq	%rdx, (%rbx)
	jmp	.L91
.L13:
	movslq	-8(%rbp), %rdx
	salq	$3, %rdx
	addq	%rdx, (%rbx)
	jmp	.L83
.L21:
	movslq	-16(%rbp), %rdx
	salq	$3, %rdx
	addq	%rdx, (%rbx)
	jmp	.L69
.L17:
	movslq	-12(%rbp), %rdx
	salq	$3, %rdx
	addq	%rdx, (%rbx)
	jmp	.L76
	.cfi_endproc
.LFE0:
	.size	getSlice, .-getSlice
	.p2align 4,,15
	.globl	getHalo
	.type	getHalo, @function
getHalo:
.LFB1:
	.cfi_startproc
	pushq	%r14
	.cfi_def_cfa_offset 16
	.cfi_offset 14, -16
	movslq	%r8d, %r8
	pushq	%r13
	.cfi_def_cfa_offset 24
	.cfi_offset 13, -24
	pushq	%r12
	.cfi_def_cfa_offset 32
	.cfi_offset 12, -32
	leaq	-8(%rdx,%r8,8), %r12
	pushq	%rbp
	.cfi_def_cfa_offset 40
	.cfi_offset 6, -40
	leaq	-4(%rcx,%r8,4), %rbp
	pushq	%rbx
	.cfi_def_cfa_offset 48
	.cfi_offset 3, -48
	subq	$16, %rsp
	.cfi_def_cfa_offset 64
	movq	%rsi, 8(%rsp)
	movq	(%r12), %rsi
	testq	%rsi, %rsi
	jle	.L105
	movq	8(%rsp), %rax
	leaq	0(,%rsi,8), %rcx
	xorl	%edx, %edx
	.p2align 4,,10
	.p2align 3
.L107:
	movsd	(%rax,%rdx), %xmm0
	movsd	%xmm0, (%rdi,%rdx)
	addq	$8, %rdx
	cmpq	%rcx, %rdx
	jne	.L107
	addq	%rdx, %rdi
	addq	%rax, %rdx
	cmpq	$1, %rsi
	movq	%rdx, 8(%rsp)
	je	.L108
.L105:
	movl	0(%rbp), %eax
	testl	%eax, %eax
	jle	.L115
	leaq	-4(%rbp), %r14
	leaq	-8(%r12), %r13
	xorl	%ebx, %ebx
	.p2align 4,,10
	.p2align 3
.L111:
	leaq	8(%rsp), %rsi
	movq	%r14, %rcx
	movq	%r13, %rdx
	addl	$1, %ebx
	call	getSlice
	cmpl	0(%rbp), %ebx
	movq	%rax, %rdi
	jl	.L111
	movq	(%r12), %rax
.L109:
	testq	%rax, %rax
	jle	.L104
.L112:
	movq	8(%rsp), %rcx
	salq	$3, %rax
	xorl	%edx, %edx
	.p2align 4,,10
	.p2align 3
.L114:
	movsd	(%rcx,%rdx), %xmm0
	movsd	%xmm0, (%rdi,%rdx)
	addq	$8, %rdx
	cmpq	%rax, %rdx
	jne	.L114
.L104:
	addq	$16, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 48
	popq	%rbx
	.cfi_def_cfa_offset 40
	popq	%rbp
	.cfi_def_cfa_offset 32
	popq	%r12
	.cfi_def_cfa_offset 24
	popq	%r13
	.cfi_def_cfa_offset 16
	popq	%r14
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L108:
	.cfi_restore_state
	movslq	0(%rbp), %rax
	salq	$3, %rax
	addq	%rax, 8(%rsp)
	movl	$1, %eax
	jmp	.L112
.L115:
	movq	%rsi, %rax
	jmp	.L109
	.cfi_endproc
.LFE1:
	.size	getHalo, .-getHalo
	.ident	"GCC: (Ubuntu 4.8.2-19ubuntu1) 4.8.2"
	.section	.note.GNU-stack,"",@progbits
