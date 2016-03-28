include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

ex1: ex1.o  chkopts
	-${CLINKER} -o ex1 ex1.o  ${PETSC_KSP_LIB}
	${RM} ex1.o


ex2: ex2.o  chkopts
	-${CLINKER} -o ex2 ex2.o  ${PETSC_KSP_LIB}
	${RM} ex2.o




#----------------------------------------------------------------------------
runex1:
	-@${MPIEXEC} -n 1 ./ex1 -ksp_monitor_short -ksp_gmres_cgs_refinement_type refine_always > ex1_1.tmp 2>&1;	  \
	   if (${DIFF} output/ex1_1.out ex1_1.tmp) then true; \
	   else printf "${PWD}\nPossible problem with ex1_1, diffs above\n=========================================\n"; fi; \
	   ${RM} -f ex1_1.tmp
runex1_2:
	-@${MPIEXEC} -n 1 ./ex1 -pc_type sor -pc_sor_symmetric -ksp_monitor_short -ksp_gmres_cgs_refinement_type refine_always >\
	   ex1_2.tmp 2>&1;   \
	   if (${DIFF} output/ex1_2.out ex1_2.tmp) then true; \
	   else printf "${PWD}\nPossible problem with ex1_2, diffs above\n=========================================\n"; fi; \
	   ${RM} -f ex1_2.tmp
runex1_3:
	-@${MPIEXEC} -n 1 ./ex1 -pc_type eisenstat -ksp_monitor_short -ksp_gmres_cgs_refinement_type refine_always >\
	   ex1_3.tmp 2>&1;   \
	   if (${DIFF} output/ex1_3.out ex1_3.tmp) then true; \
	   else printf "${PWD}\nPossible problem with ex1_3, diffs above\n=========================================\n"; fi; \

runex2:
	-@${MPIEXEC} -n 1 ./ex2 -ksp_monitor_short -m 5 -n 5 -ksp_gmres_cgs_refinement_type refine_always > ex2_1.tmp 2>&1; \
	   if (${DIFF} output/ex2_1.out ex2_1.tmp) then true; \
	   else printf "${PWD}\nPossible problem with ex2_1, diffs above\n=========================================\n"; fi; \
	   ${RM} -f ex2_1.tmp
runex2_2:
	-@${MPIEXEC} -n 2 ./ex2 -ksp_monitor_short -m 5 -n 5 -ksp_gmres_cgs_refinement_type refine_always > ex2_2.tmp 2>&1; \
	   if (${DIFF} output/ex2_2.out ex2_2.tmp) then true; \
	   else printf "${PWD}\nPossible problem with ex2_2, diffs above\n=========================================\n"; fi; \
	   ${RM} -f ex2_2.tmp
runex2_3:
	-@${MPIEXEC} -n 1 ./ex2 -pc_type sor -pc_sor_symmetric -ksp_monitor_short -ksp_gmres_cgs_refinement_type refine_always > \
	    ex2_3.tmp 2>&1;   \
	   if (${DIFF} output/ex2_3.out ex2_3.tmp) then true; \
	   else printf "${PWD}\nPossible problem with ex2_3, diffs above\n=========================================\n"; fi; \
	   ${RM} -f ex2_3.tmp
runex2_4:
	-@${MPIEXEC} -n 1 ./ex2 -pc_type eisenstat -ksp_monitor_short -ksp_gmres_cgs_refinement_type refine_always >\
	    ex2_4.tmp 2>&1;   \
	   if (${DIFF} output/ex2_4.out ex2_4.tmp) then true; \
	   else printf "${PWD}\nPossible problem with ex2_4, diffs above\n=========================================\n"; fi; \
	   ${RM} -f ex2_4.tmp
runex2_5:
	-@${MPIEXEC} -n 2 ./ex2 -ksp_monitor_short -m 5 -n 5 -mat_view draw -ksp_gmres_cgs_refinement_type refine_always -nox  > ex2_5.tmp 2>&1; \
	   if (${DIFF} output/ex2_2.out ex2_5.tmp) then true; \
	   else printf "${PWD}\nPossible problem with ex2_5, diffs above\n=========================================\n"; fi; \
	   ${RM} -f ex2_5.tmp
runex2_bjacobi:
	-@${MPIEXEC} -n 4 ./ex2 -pc_type bjacobi -pc_bjacobi_blocks 1 -ksp_monitor_short -sub_pc_type jacobi -sub_ksp_type gmres > ex2.tmp 2>&1; \
	   if (${DIFF} output/ex2_bjacobi.out ex2.tmp) then true; \
	   else printf "${PWD}\nPossible problem with ex2_bjacobi, diffs above\n=========================================\n"; fi; \
	   ${RM} -f ex2.tmp
runex2_bjacobi_2:
	-@${MPIEXEC} -n 4 ./ex2 -pc_type bjacobi -pc_bjacobi_blocks 2 -ksp_monitor_short -sub_pc_type jacobi -sub_ksp_type gmres -ksp_view > ex2.tmp 2>&1; \
	   if (${DIFF} output/ex2_bjacobi_2.out ex2.tmp) then true; \
	   else printf "${PWD}\nPossible problem with ex2_bjacobi_2, diffs above\n=========================================\n"; fi; \
	   ${RM} -f ex2.tmp
runex2_bjacobi_3:
	-@${MPIEXEC} -n 4 ./ex2 -pc_type bjacobi -pc_bjacobi_blocks 4 -ksp_monitor_short -sub_pc_type jacobi -sub_ksp_type gmres > ex2.tmp 2>&1; \
	   if (${DIFF} output/ex2_bjacobi_3.out ex2.tmp) then true; \
	   else printf "${PWD}\nPossible problem with ex2_bjacobi_3, diffs above\n=========================================\n"; fi; \
	   ${RM} -f ex2.tmp
runex2_chebyest_1:
	-@${MPIEXEC} -n 1 ./ex2 -m 80 -n 80 -ksp_pc_side right -pc_type ksp -ksp_ksp_type chebyshev -ksp_ksp_max_it 5 -ksp_ksp_chebyshev_esteig 0.9,0,0,1.1 -ksp_monitor_short > ex2.tmp 2>&1; \
           ${DIFF} output/ex2_chebyest_1.out ex2.tmp || printf "${PWD}\nPossible problem with ex2_chebyest_1, diffs above\n=========================================\n"; \
           ${RM} -f ex2.tmp
runex2_chebyest_2:
	-@${MPIEXEC} -n 1 ./ex2 -m 80 -n 80 -ksp_pc_side right -pc_type ksp -ksp_ksp_type chebyshev -ksp_ksp_max_it 5 -ksp_ksp_chebyshev_esteig 0.9,0,0,1.1 -ksp_esteig_ksp_type cg -ksp_monitor_short > ex2.tmp 2>&1; \
           ${DIFF} output/ex2_chebyest_2.out ex2.tmp || printf "${PWD}\nPossible problem with ex2_chebyest_2, diffs above\n=========================================\n"; \
           ${RM} -f ex2.tmp
runex2_umfpack:
	-@${MPIEXEC} -n 1 ./ex2 -ksp_type preonly -pc_type lu -pc_factor_mat_solver_package umfpack > ex2_umfpack.tmp 2>&1; \
           if (${DIFF} output/ex2_umfpack.out ex2_umfpack.tmp) then true; \
           else printf "${PWD}\nPossible problem with ex2_umfpack, diffs above\n=========================================\n"; fi; \
           ${RM} -f ex2_umfpack.tmp
runex2_mkl_pardiso_lu:
	-@${MPIEXEC} -n 1 ./ex2 -ksp_type preonly -pc_type lu -pc_factor_mat_solver_package mkl_pardiso > ex2_mkl_pardiso.tmp 2>&1; \
           if (${DIFF} output/ex2_mkl_pardiso_lu.out ex2_mkl_pardiso.tmp) then true; \
           else printf "${PWD}\nPossible problem with ex2_mkl_pardiso_lu, diffs above\n=========================================\n"; fi; \
           ${RM} -f ex2_mkl_pardiso.tmp
runex2_mkl_pardiso_cholesky:
	-@${MPIEXEC} -n 1 ./ex2 -ksp_type preonly -pc_type cholesky -mat_type sbaij -pc_factor_mat_solver_package mkl_pardiso > ex2_mkl_pardiso.tmp 2>&1; \
           if (${DIFF} output/ex2_mkl_pardiso_cholesky.out ex2_mkl_pardiso.tmp) then true; \
           else printf "${PWD}\nPossible problem with ex2_mkl_pardiso_cholesky, diffs above\n=========================================\n"; fi; \
           ${RM} -f ex2_mkl_pardiso.tmp
runex2_fbcgs:
	-@${MPIEXEC} -n 1 ./ex2 -ksp_type fbcgs -pc_type ilu  > ex2.tmp 2>&1; \
           if (${DIFF} output/ex2_fbcgs.out ex2.tmp) then true; \
           else printf "${PWD}\nPossible problem with ex2_fbcgs, diffs above\n=========================================\n"; fi; \
           ${RM} -f ex2.tmp
runex2_fbcgs_2:
	-@${MPIEXEC} -n 3 ./ex2 -ksp_type fbcgsr -pc_type bjacobi > ex2.tmp 2>&1; \
           if (${DIFF} output/ex2_fbcgs_2.out ex2.tmp) then true; \
           else printf "${PWD}\nPossible problem with ex2_fbcgs_2, diffs above\n=========================================\n"; fi; \
           ${RM} -f ex2.tmp
runex2_pipecg:
	-@${MPIEXEC} -n 1 ./ex2 -ksp_monitor_short -ksp_type pipecg -m 9 -n 9 > ex2_pipecg.tmp 2>&1; \
	   ${DIFF} output/ex2_pipecg.out ex2_pipecg.tmp || printf "${PWD}\nPossible problem with ex2_pipecg, diffs above\n=========================================\n"; \
	   ${RM} -f ex2_pipecg.tmp
runex2_pipecr:
	-@${MPIEXEC} -n 1 ./ex2 -ksp_monitor_short -ksp_type pipecr -m 9 -n 9 > ex2_pipecr.tmp 2>&1; \
	   ${DIFF} output/ex2_pipecr.out ex2_pipecr.tmp || printf "${PWD}\nPossible problem with ex2_pipecr, diffs above\n=========================================\n"; \
	   ${RM} -f ex2_pipecr.tmp
runex2_groppcg:
	-@${MPIEXEC} -n 1 ./ex2 -ksp_monitor_short -ksp_type groppcg -m 9 -n 9 > ex2_groppcg.tmp 2>&1; \
	   ${DIFF} output/ex2_groppcg.out ex2_groppcg.tmp || printf "${PWD}\nPossible problem with ex2_groppcg, diffs above\n=========================================\n"; \
	   ${RM} -f ex2_groppcg.tmp

