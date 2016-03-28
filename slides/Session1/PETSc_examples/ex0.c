
static char help[] = "Solves a linear system in sequential with KSP and LU decomposition as preconditioner.\n\n";
/*T
Concepts: KSP^solving a system of linear equations
Processors: 1
T*/

/*
   Include "petscksp.h" so that we can use KSP solvers.  Note that this file
   automatically includes:
   petscsys.h       - base PETSc routines   petscvec.h - vectors
   petscmat.h - matrices
   petscis.h     - index sets            petscksp.h - Krylov subspace methods
   petscviewer.h - viewers               petscpc.h  - preconditioners

Note:  The corresponding parallel example is ex23.c
*/
#include <petscksp.h>

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{
	Vec            x, b, u;      /* approx solution, RHS, exact solution */
	Mat            A;            /* linear system matrix */
	KSP            ksp;         /* linear solver context */
	PC             pc;           /* preconditioner context */
	PetscReal      norm,tol=1.e-14;  /* norm of solution error */
	PetscErrorCode ierr;
	PetscInt       n = 10,its;
	PetscMPIInt    size;
	PetscScalar    neg_one      = -1.0,one = 1.0;
	PetscBool      nonzeroguess = PETSC_FALSE;


	PetscRandom    rctx;     /* random number generator context */

	PetscInitialize(&argc,&args,(char*)0,help);
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
	if (size != 1) SETERRQ(PETSC_COMM_WORLD,1,"This is a uniprocessor example only!");
	ierr = PetscOptionsGetInt(NULL,"-n",&n,NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetBool(NULL,"-nonzero_guess",&nonzeroguess,NULL);CHKERRQ(ierr);


	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	   Compute the matrix and right-hand-side vector that define
	   the linear system, Ax = b.
	   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

	/*
	   Create vectors.  Note that we form 1 vector from scratch and
	   then duplicate as needed.
	   */
	ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
	ierr = PetscObjectSetName((PetscObject) x, "Solution");CHKERRQ(ierr);
	ierr = VecSetSizes(x,PETSC_DECIDE,n);CHKERRQ(ierr);
	ierr = VecSetFromOptions(x);CHKERRQ(ierr);
	ierr = VecDuplicate(x,&b);CHKERRQ(ierr);
	ierr = VecDuplicate(x,&u);CHKERRQ(ierr);

	/*
	   Create matrix. 
	   */
	ierr = MatCreateSeqDense(PETSC_COMM_SELF,n,n,PETSC_NULL,&A);CHKERRQ(ierr);

	ierr = PetscRandomCreate(PETSC_COMM_WORLD,&rctx);CHKERRQ(ierr);
	ierr = PetscRandomSetFromOptions(rctx);CHKERRQ(ierr);
	ierr =MatSetRandom(A,rctx);CHKERRQ(ierr);
	ierr =PetscRandomDestroy(&rctx);CHKERRQ(ierr);



	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

	//MatView(A,PETSC_VIEWER_STDOUT_WORLD);

	/*
	   Set exact solution; then compute right-hand-side vector.
	   */
	ierr = VecSet(u,one);CHKERRQ(ierr);
	ierr = MatMult(A,u,b);CHKERRQ(ierr);

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	   Create the linear solver and set various options
	   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	/*
	   Create linear solver context
	   */
	ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);

	/*
	   Set operators. Here the matrix that defines the linear system
	   also serves as the preconditioning matrix.
	   */
	ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);

	/*
	   Set linear solver defaults for this problem (optional).
	   - By extracting the KSP and PC contexts from the KSP context,
	   we can then directly call any KSP and PC routines to set
	   various options.
	   - The following four statements are optional; all of these
	   parameters could alternatively be specified at runtime via
	   KSPSetFromOptions();
	   */
	ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
	ierr = PCSetType(pc,PCJACOBI);CHKERRQ(ierr);
	ierr = KSPSetTolerances(ksp,1.e-5,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);

	/*
	   Set runtime options, e.g.,
	   -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
	   These options will override those specified above as long as
	   KSPSetFromOptions() is called _after_ any other customization
	   routines.
	   */
	ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

	if (nonzeroguess) {
		PetscScalar p = .5;
		ierr = VecSet(x,p);CHKERRQ(ierr);
		ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);CHKERRQ(ierr);
	}
    ierr=KSPSetUp(ksp); CHKERRQ(ierr);
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	   Solve the linear system
	   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	/*
	   Solve linear system
	   */
	ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);

	/*
	   View solver info; we could instead use the option -ksp_view to
	   print this info to the screen at the conclusion of KSPSolve().
	   */
	ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	   Check solution and clean up
	   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	/*
	   Check the error
	   */
	ierr = VecAXPY(x,neg_one,u);CHKERRQ(ierr);
	ierr = VecNorm(x,NORM_2,&norm);CHKERRQ(ierr);
	ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
	if (norm > tol) {
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of error %g, Iterations %D\n",(double)norm,its);CHKERRQ(ierr);
	}

	/*
	   Free work space.  All PETSc objects should be destroyed when they
	   are no longer needed.
	   */
	ierr = VecDestroy(&x);CHKERRQ(ierr); ierr = VecDestroy(&u);CHKERRQ(ierr);
	ierr = VecDestroy(&b);CHKERRQ(ierr); ierr = MatDestroy(&A);CHKERRQ(ierr);
	ierr = KSPDestroy(&ksp);CHKERRQ(ierr);

	/*
	   Always call PetscFinalize() before exiting a program.  This routine
	   - finalizes the PETSc libraries as well as MPI
	   - provides summary and diagnostic information if certain runtime
	   options are chosen (e.g., -log_summary).
	   */
	ierr = PetscFinalize();
	return 0;
}
