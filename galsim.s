	.file	"galsim.c"
	.section	.text.unlikely,"ax",@progbits
.LCOLDB3:
	.text
.LHOTB3:
	.p2align 4,,15
	.globl	force
	.type	force, @function
force:
.LFB53:
	.cfi_startproc
	testl	%esi, %esi
	jle	.L19
	pushq	%r15
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
	pushq	%r14
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	leal	-1(%rsi), %edx
	pushq	%r13
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	pushq	%r12
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	leaq	80(%rdi), %r12
	pushq	%rbp
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	pushq	%rbx
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	movq	%rdi, %r14
	pxor	%xmm8, %xmm8
	movl	$1, %r15d
	subq	$72, %rsp
	.cfi_def_cfa_offset 128
	xorl	%r13d, %r13d
	movsd	.LC1(%rip), %xmm7
	movsd	.LC2(%rip), %xmm5
	xorpd	%xmm0, %xmm7
	movapd	%xmm5, %xmm6
	.p2align 4,,10
	.p2align 3
.L12:
	testl	%r13d, %r13d
	movsd	-64(%r12), %xmm3
	movsd	-56(%r12), %xmm4
	jle	.L13
	leal	-1(%r13), %eax
	leaq	16(%r14), %rbx
	movapd	%xmm8, %xmm2
	salq	$6, %rax
	leaq	80(%r14,%rax), %rbp
	.p2align 4,,10
	.p2align 3
.L6:
	movapd	%xmm3, %xmm9
	movapd	%xmm4, %xmm1
	subsd	(%rbx), %xmm9
	subsd	8(%rbx), %xmm1
	movapd	%xmm9, %xmm0
	mulsd	%xmm1, %xmm1
	mulsd	%xmm9, %xmm0
	addsd	%xmm1, %xmm0
	sqrtsd	%xmm0, %xmm1
	ucomisd	%xmm1, %xmm1
	jp	.L22
.L4:
	addsd	%xmm5, %xmm1
	addq	$64, %rbx
	movapd	%xmm1, %xmm0
	mulsd	%xmm1, %xmm0
	mulsd	%xmm0, %xmm1
	movsd	-32(%rbx), %xmm0
	cmpq	%rbp, %rbx
	divsd	%xmm1, %xmm0
	mulsd	%xmm0, %xmm9
	addsd	%xmm9, %xmm2
	jne	.L6
.L3:
	addl	$1, %r13d
	cmpl	%r13d, %esi
	je	.L7
	movl	%edx, %eax
	movq	%r12, %rbx
	subl	%r13d, %eax
	addq	%r15, %rax
	salq	$6, %rax
	leaq	80(%r14,%rax), %rbp
	.p2align 4,,10
	.p2align 3
.L10:
	movapd	%xmm3, %xmm9
	movapd	%xmm4, %xmm1
	subsd	(%rbx), %xmm9
	subsd	8(%rbx), %xmm1
	movapd	%xmm9, %xmm0
	mulsd	%xmm1, %xmm1
	mulsd	%xmm9, %xmm0
	addsd	%xmm1, %xmm0
	sqrtsd	%xmm0, %xmm1
	ucomisd	%xmm1, %xmm1
	jp	.L23
.L8:
	addsd	%xmm6, %xmm1
	addq	$64, %rbx
	movapd	%xmm1, %xmm0
	mulsd	%xmm1, %xmm0
	mulsd	%xmm0, %xmm1
	movsd	-32(%rbx), %xmm0
	cmpq	%rbx, %rbp
	divsd	%xmm1, %xmm0
	mulsd	%xmm0, %xmm9
	addsd	%xmm9, %xmm2
	jne	.L10
	movsd	-32(%r12), %xmm0
	addq	$1, %r15
	addq	$64, %r12
	mulsd	%xmm7, %xmm0
	mulsd	%xmm0, %xmm2
	movsd	%xmm2, -112(%r12)
	movsd	%xmm2, -104(%r12)
	jmp	.L12
.L23:
	movl	%esi, 60(%rsp)
	movsd	%xmm6, 48(%rsp)
	movl	%edx, 56(%rsp)
	movsd	%xmm5, 40(%rsp)
	movsd	%xmm8, 32(%rsp)
	movsd	%xmm7, 24(%rsp)
	movsd	%xmm2, 16(%rsp)
	movsd	%xmm4, 8(%rsp)
	movsd	%xmm3, (%rsp)
	call	sqrt
	movsd	(%rsp), %xmm3
	movapd	%xmm0, %xmm1
	movl	60(%rsp), %esi
	movl	56(%rsp), %edx
	movapd	%xmm3, %xmm9
	movsd	48(%rsp), %xmm6
	subsd	(%rbx), %xmm9
	movsd	40(%rsp), %xmm5
	movsd	32(%rsp), %xmm8
	movsd	24(%rsp), %xmm7
	movsd	16(%rsp), %xmm2
	movsd	8(%rsp), %xmm4
	jmp	.L8
.L22:
	movl	%esi, 60(%rsp)
	movsd	%xmm6, 48(%rsp)
	movl	%edx, 56(%rsp)
	movsd	%xmm5, 40(%rsp)
	movsd	%xmm8, 32(%rsp)
	movsd	%xmm7, 24(%rsp)
	movsd	%xmm2, 16(%rsp)
	movsd	%xmm4, 8(%rsp)
	movsd	%xmm3, (%rsp)
	call	sqrt
	movsd	(%rsp), %xmm3
	movapd	%xmm0, %xmm1
	movl	60(%rsp), %esi
	movl	56(%rsp), %edx
	movapd	%xmm3, %xmm9
	movsd	48(%rsp), %xmm6
	subsd	(%rbx), %xmm9
	movsd	40(%rsp), %xmm5
	movsd	32(%rsp), %xmm8
	movsd	24(%rsp), %xmm7
	movsd	16(%rsp), %xmm2
	movsd	8(%rsp), %xmm4
	jmp	.L4
.L7:
	salq	$6, %rdx
	addq	%rdx, %r14
	mulsd	48(%r14), %xmm7
	mulsd	%xmm7, %xmm2
	movsd	%xmm2, 32(%r14)
	movsd	%xmm2, 40(%r14)
	addq	$72, %rsp
	.cfi_def_cfa_offset 56
	popq	%rbx
	.cfi_restore 3
	.cfi_def_cfa_offset 48
	popq	%rbp
	.cfi_restore 6
	.cfi_def_cfa_offset 40
	popq	%r12
	.cfi_restore 12
	.cfi_def_cfa_offset 32
	popq	%r13
	.cfi_restore 13
	.cfi_def_cfa_offset 24
	popq	%r14
	.cfi_restore 14
	.cfi_def_cfa_offset 16
	popq	%r15
	.cfi_restore 15
	.cfi_def_cfa_offset 8
.L19:
	rep ret
.L13:
	.cfi_def_cfa_offset 128
	.cfi_offset 3, -56
	.cfi_offset 6, -48
	.cfi_offset 12, -40
	.cfi_offset 13, -32
	.cfi_offset 14, -24
	.cfi_offset 15, -16
	movapd	%xmm8, %xmm2
	jmp	.L3
	.cfi_endproc
.LFE53:
	.size	force, .-force
	.section	.text.unlikely
.LCOLDE3:
	.text
.LHOTE3:
	.section	.text.unlikely
.LCOLDB4:
	.text
.LHOTB4:
	.p2align 4,,15
	.globl	timestep
	.type	timestep, @function
timestep:
.LFB54:
	.cfi_startproc
	testl	%esi, %esi
	jle	.L24
	leal	-1(%rsi), %eax
	addq	$1, %rax
	salq	$6, %rax
	addq	%rdi, %rax
	.p2align 4,,10
	.p2align 3
.L26:
	movsd	48(%rdi), %xmm3
	addq	$64, %rdi
	movsd	-24(%rdi), %xmm1
	movsd	-32(%rdi), %xmm2
	divsd	%xmm3, %xmm1
	divsd	%xmm3, %xmm2
	mulsd	%xmm0, %xmm1
	mulsd	%xmm0, %xmm2
	addsd	-56(%rdi), %xmm1
	addsd	-64(%rdi), %xmm2
	movsd	%xmm1, -56(%rdi)
	mulsd	%xmm0, %xmm1
	movsd	%xmm2, -64(%rdi)
	mulsd	%xmm0, %xmm2
	addsd	-40(%rdi), %xmm1
	addsd	-48(%rdi), %xmm2
	movsd	%xmm1, -40(%rdi)
	movsd	%xmm2, -48(%rdi)
	cmpq	%rdi, %rax
	jne	.L26
.L24:
	rep ret
	.cfi_endproc
.LFE54:
	.size	timestep, .-timestep
	.section	.text.unlikely
.LCOLDE4:
	.text
.LHOTE4:
	.section	.rodata.str1.8,"aMS",@progbits,1
	.align 8
.LC6:
	.string	"ERROR! Expected five input arguments: N filename nsteps delta_t graphics"
	.section	.rodata.str1.1,"aMS",@progbits,1
.LC7:
	.string	"rb"
.LC8:
	.string	"FATAL ERROR!"
.LC11:
	.string	"wb"
.LC12:
	.string	"result.gal"
	.section	.rodata.str1.8
	.align 8
.LC14:
	.string	"computation time: %7.3f wall seconds\n"
	.section	.text.unlikely
.LCOLDB15:
	.section	.text.startup,"ax",@progbits
.LHOTB15:
	.p2align 4,,15
	.globl	main
	.type	main, @function
main:
.LFB55:
	.cfi_startproc
	pushq	%r15
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
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
	movl	%edi, %ebx
	movq	%rsi, %rbp
	xorl	%esi, %esi
	subq	$104, %rsp
	.cfi_def_cfa_offset 160
	leaq	64(%rsp), %rdi
	movq	%fs:40, %rax
	movq	%rax, 88(%rsp)
	xorl	%eax, %eax
	call	gettimeofday
	pxor	%xmm0, %xmm0
	cmpl	$6, %ebx
	pxor	%xmm1, %xmm1
	cvtsi2sdq	72(%rsp), %xmm0
	divsd	.LC5(%rip), %xmm0
	cvtsi2sdq	64(%rsp), %xmm1
	addsd	%xmm1, %xmm0
	movsd	%xmm0, 40(%rsp)
	jne	.L66
	movq	8(%rbp), %rdi
	xorl	%esi, %esi
	movl	$10, %edx
	call	strtol
	movq	24(%rbp), %rdi
	movl	$10, %edx
	xorl	%esi, %esi
	movq	%rax, %r14
	movq	16(%rbp), %rbx
	movq	%rax, 32(%rsp)
	movl	%eax, %r13d
	call	strtol
	movq	32(%rbp), %rdi
	xorl	%esi, %esi
	movq	%rax, 56(%rsp)
	movl	%eax, 24(%rsp)
	call	strtod
	movq	40(%rbp), %rdi
	xorl	%esi, %esi
	movl	$10, %edx
	movsd	%xmm0, 8(%rsp)
	call	strtol
	movq	%rax, 48(%rsp)
	movl	%eax, 28(%rsp)
	movl	$100, %eax
	cltd
	pxor	%xmm7, %xmm7
	idivl	%r14d
	movslq	%r14d, %rdi
	salq	$6, %rdi
	cvtsi2sd	%eax, %xmm7
	movsd	%xmm7, 16(%rsp)
	call	malloc
	movl	$.LC7, %esi
	movq	%rbx, %rdi
	movq	%rax, %r12
	call	fopen
	testl	%r14d, %r14d
	movq	%rax, %rbp
	movsd	8(%rsp), %xmm6
	jle	.L36
	movq	32(%rsp), %rax
	leaq	16(%r12), %rbx
	subl	$1, %eax
	salq	$6, %rax
	leaq	80(%r12,%rax), %r15
	.p2align 4,,10
	.p2align 3
.L35:
	movq	%rbx, %rdi
	movq	%rbp, %rcx
	movl	$1, %edx
	movl	$8, %esi
	movsd	%xmm6, 8(%rsp)
	call	fread
	leaq	8(%rbx), %rdi
	movq	%rbp, %rcx
	movl	$1, %edx
	movl	$8, %esi
	call	fread
	leaq	32(%rbx), %rdi
	movq	%rbp, %rcx
	movl	$1, %edx
	movl	$8, %esi
	call	fread
	leaq	-16(%rbx), %rdi
	movq	%rbp, %rcx
	movl	$1, %edx
	movl	$8, %esi
	call	fread
	leaq	-8(%rbx), %rdi
	movq	%rbp, %rcx
	movl	$1, %edx
	movl	$8, %esi
	call	fread
	leaq	40(%rbx), %rdi
	movq	%rbp, %rcx
	movl	$1, %edx
	movl	$8, %esi
	addq	$64, %rbx
	call	fread
	cmpq	%r15, %rbx
	movsd	8(%rsp), %xmm6
	jne	.L35
.L36:
	cmpl	$1, 48(%rsp)
	je	.L67
.L34:
	movl	56(%rsp), %edx
	testl	%edx, %edx
	jle	.L45
	movq	32(%rsp), %rax
	xorl	%ebx, %ebx
	subl	$1, %eax
	salq	$6, %rax
	leaq	72(%r12,%rax), %r15
	leaq	88(%r12,%rax), %rbp
	.p2align 4,,10
	.p2align 3
.L44:
	movsd	16(%rsp), %xmm0
	movl	%r13d, %esi
	movq	%r12, %rdi
	movsd	%xmm6, 8(%rsp)
	call	force
	testl	%r13d, %r13d
	movsd	8(%rsp), %xmm6
	jle	.L40
	leaq	8(%r12), %rax
	.p2align 4,,10
	.p2align 3
.L41:
	movsd	40(%rax), %xmm2
	addq	$64, %rax
	movsd	-32(%rax), %xmm0
	movsd	-40(%rax), %xmm1
	divsd	%xmm2, %xmm0
	divsd	%xmm2, %xmm1
	mulsd	%xmm6, %xmm0
	mulsd	%xmm6, %xmm1
	addsd	-64(%rax), %xmm0
	addsd	-72(%rax), %xmm1
	movsd	%xmm0, -64(%rax)
	mulsd	%xmm6, %xmm0
	movsd	%xmm1, -72(%rax)
	mulsd	%xmm6, %xmm1
	addsd	-48(%rax), %xmm0
	addsd	-56(%rax), %xmm1
	movsd	%xmm0, -48(%rax)
	movsd	%xmm1, -56(%rax)
	cmpq	%rax, %r15
	jne	.L41
	cmpl	$1, 28(%rsp)
	je	.L68
.L49:
	addl	$1, %ebx
	cmpl	24(%rsp), %ebx
	jne	.L44
.L45:
	movl	$.LC11, %esi
	movl	$.LC12, %edi
	call	fopen
	movq	%rax, %rbp
	movl	32(%rsp), %eax
	testl	%eax, %eax
	jle	.L39
	movq	32(%rsp), %rax
	leaq	16(%r12), %rbx
	subl	$1, %eax
	salq	$6, %rax
	leaq	80(%r12,%rax), %r12
	.p2align 4,,10
	.p2align 3
.L46:
	movq	%rbx, %rdi
	movq	%rbp, %rcx
	movl	$1, %edx
	movl	$8, %esi
	call	fwrite
	leaq	8(%rbx), %rdi
	movq	%rbp, %rcx
	movl	$1, %edx
	movl	$8, %esi
	call	fwrite
	leaq	32(%rbx), %rdi
	movq	%rbp, %rcx
	movl	$1, %edx
	movl	$8, %esi
	call	fwrite
	leaq	-16(%rbx), %rdi
	movq	%rbp, %rcx
	movl	$1, %edx
	movl	$8, %esi
	call	fwrite
	leaq	-8(%rbx), %rdi
	movq	%rbp, %rcx
	movl	$1, %edx
	movl	$8, %esi
	call	fwrite
	leaq	40(%rbx), %rdi
	movq	%rbp, %rcx
	movl	$1, %edx
	movl	$8, %esi
	addq	$64, %rbx
	call	fwrite
	cmpq	%r12, %rbx
	jne	.L46
.L39:
	leaq	64(%rsp), %rdi
	xorl	%esi, %esi
	call	gettimeofday
	pxor	%xmm0, %xmm0
	movl	$.LC14, %esi
	pxor	%xmm1, %xmm1
	movl	$1, %edi
	movl	$1, %eax
	cvtsi2sdq	72(%rsp), %xmm0
	divsd	.LC5(%rip), %xmm0
	cvtsi2sdq	64(%rsp), %xmm1
	addsd	%xmm1, %xmm0
	subsd	40(%rsp), %xmm0
	call	__printf_chk
	xorl	%eax, %eax
.L31:
	movq	88(%rsp), %rcx
	xorq	%fs:40, %rcx
	jne	.L69
	addq	$104, %rsp
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
.L68:
	.cfi_restore_state
	movsd	%xmm6, 8(%rsp)
	leaq	24(%r12), %r14
	call	ClearScreen
	movsd	8(%rsp), %xmm6
	.p2align 4,,10
	.p2align 3
.L43:
	pxor	%xmm0, %xmm0
	movss	.LC9(%rip), %xmm3
	pxor	%xmm1, %xmm1
	addq	$64, %r14
	pxor	%xmm5, %xmm5
	movsd	%xmm6, 8(%rsp)
	movaps	%xmm3, %xmm2
	cvtsd2ss	-72(%r14), %xmm0
	cvtsd2ss	-64(%r14), %xmm1
	movss	.LC13(%rip), %xmm4
	call	DrawCircle
	cmpq	%rbp, %r14
	movsd	8(%rsp), %xmm6
	jne	.L43
.L47:
	movsd	%xmm6, 8(%rsp)
	call	Refresh
	movl	$3000, %edi
	call	usleep
	movsd	8(%rsp), %xmm6
	jmp	.L49
.L40:
	cmpl	$1, 28(%rsp)
	jne	.L49
	movsd	%xmm6, 8(%rsp)
	call	ClearScreen
	movsd	8(%rsp), %xmm6
	jmp	.L47
.L67:
	movl	$800, %edx
	movl	$1400, %esi
	movl	$.LC8, %edi
	movsd	%xmm6, 8(%rsp)
	call	InitializeGraphics
	pxor	%xmm0, %xmm0
	movss	.LC9(%rip), %xmm1
	call	SetCAxes
	movsd	8(%rsp), %xmm6
	jmp	.L34
.L66:
	movl	$.LC6, %edi
	call	puts
	orl	$-1, %eax
	jmp	.L31
.L69:
	call	__stack_chk_fail
	.cfi_endproc
.LFE55:
	.size	main, .-main
	.section	.text.unlikely
.LCOLDE15:
	.section	.text.startup
.LHOTE15:
	.section	.rodata.cst16,"aM",@progbits,16
	.align 16
.LC1:
	.long	0
	.long	-2147483648
	.long	0
	.long	0
	.section	.rodata.cst8,"aM",@progbits,8
	.align 8
.LC2:
	.long	3539053052
	.long	1062232653
	.align 8
.LC5:
	.long	0
	.long	1093567616
	.section	.rodata.cst4,"aM",@progbits,4
	.align 4
.LC9:
	.long	1065353216
	.align 4
.LC13:
	.long	998445679
	.ident	"GCC: (Ubuntu 5.4.0-6ubuntu1~16.04.4) 5.4.0 20160609"
	.section	.note.GNU-stack,"",@progbits
