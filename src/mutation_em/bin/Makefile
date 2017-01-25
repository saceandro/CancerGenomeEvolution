CXX = g++

generate: generate.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_wf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ generate.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

generate_params: generate_params.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_wf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ generate_params.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

generate_reads: generate_reads.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_wf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ generate_reads.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

compare: compare.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_mutation_em.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ compare.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

compare_wo_resp: compare_wo_resp.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_mutation_em.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ compare_wo_resp.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

compare_xi: compare_xi.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_mutation_em.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ compare_xi.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

compare_xi_corrected: compare_xi_corrected.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_mutation_em.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ compare_xi_corrected.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

compare_xi_log: compare_xi_log.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_mutation_em.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ compare_xi_log.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

compare_precalc_resp: compare_precalc_resp.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_mutation_em.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -Ofast -o $@ compare_precalc_resp.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

compare_precalc_resp_u: compare_precalc_resp_u.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_mutation_em.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -Ofast -o $@ compare_precalc_resp_u.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

compare_precalc_resp_n: compare_precalc_resp_n.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_mutation_em.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -Ofast -o $@ compare_precalc_resp_n.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

compare_precalc_resp_with_normal: compare_precalc_resp_with_normal.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_mutation_em.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -Ofast -o $@ compare_precalc_resp_with_normal.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

compare_precalc_resp_with_normal_u: compare_precalc_resp_with_normal_u.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_mutation_em.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -Ofast -o $@ compare_precalc_resp_with_normal_u.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

compare_precalc_resp_with_normal_n: compare_precalc_resp_with_normal_n.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_mutation_em.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -Ofast -o $@ compare_precalc_resp_with_normal_n.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

compare_omega: compare_omega.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_mutation_em.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ compare_omega.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

compare_vf_t: compare_vf_t.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_mutation_em.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ compare_vf_t.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

compare_simple: compare_simple.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_wf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ compare_simple.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

grad: grad.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_wf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ grad.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

grad_desc: grad_desc.cc ../variant_fraction.cc ../setting.hh ../../../util/loglib_header.hh ../../../util/enumtree_mutation_em_vf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -Ofast -o $@ grad_desc.cc ../variant_fraction.cc ../../../util/loglib_header.cc -lgsl -lgslcblas -lm

grad_desc_multinomial: grad_desc_multinomial.cc ../variant_fraction.cc ../setting.hh ../../../util/loglib_header.hh ../../../util/enumtree_mutation_em_multinomial.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -Ofast -o $@ grad_desc_multinomial.cc ../variant_fraction.cc ../../../util/loglib_header.cc -lgsl -lgslcblas -lm

grad_desc_llik: grad_desc_llik.cc ../variant_fraction.cc ../setting.hh ../../../util/loglib_header.hh ../../../util/enumtree_mutation_em_vf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -Ofast -o $@ grad_desc_llik.cc ../variant_fraction.cc ../../../util/loglib_header.cc -lgsl -lgslcblas -lm

grad_desc_llik_multinomial: grad_desc_llik_multinomial.cc ../variant_fraction.cc ../setting.hh ../../../util/loglib_header.hh ../../../util/enumtree_mutation_em_multinomial.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -Ofast -o $@ grad_desc_llik_multinomial.cc ../variant_fraction.cc ../../../util/loglib_header.cc -lgsl -lgslcblas -lm

grad_desc_llik_calc_multinomial: grad_desc_llik_calc_multinomial.cc ../variant_fraction.cc ../setting.hh ../../../util/loglib_header.hh ../../../util/enumtree_mutation_em_multinomial.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -Ofast -o $@ grad_desc_llik_calc_multinomial.cc ../variant_fraction.cc ../../../util/loglib_header.cc -lgsl -lgslcblas -lm

grad_foraccuracy: grad_foraccuracy.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_wf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -Ofast -o $@ grad_foraccuracy.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

grad_foraccuracy_n: grad_foraccuracy_n.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_wf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -Ofast -o $@ grad_foraccuracy_n.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

grad_foraccuracy_n_init: grad_foraccuracy_n_init.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_wf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -Ofast -o $@ grad_foraccuracy_n_init.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

em: em.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_mutation_em.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -Ofast -o $@ em.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

estep: estep.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_mutation_em_vf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -Ofast -o $@ estep.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

estep_mem: estep_mem.cc ../variant_fraction.cc ../setting.hh ../../../util/loglib_header.hh ../../../util/enumtree_mutation_em_vf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -Ofast -o $@ estep_mem.cc ../variant_fraction.cc ../../../util/loglib_header.cc -lgsl -lgslcblas -lm

estep_mem_multinomial: estep_mem_multinomial.cc ../variant_fraction.cc ../setting.hh ../../../util/loglib_header.hh ../../../util/enumtree_mutation_em_multinomial.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -Ofast -o $@ estep_mem_multinomial.cc ../variant_fraction.cc ../../../util/loglib_header.cc -lgsl -lgslcblas -lm

estep_mem_llik: estep_mem_llik.cc ../variant_fraction.cc ../setting.hh ../../../util/loglib_header.hh ../../../util/enumtree_mutation_em_vf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -Ofast -o $@ estep_mem_llik.cc ../variant_fraction.cc ../../../util/loglib_header.cc -lgsl -lgslcblas -lm

estep_mem_resp: estep_mem_resp.cc ../variant_fraction.cc ../setting.hh ../../../util/loglib_header.hh ../../../util/enumtree_mutation_em_vf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -Ofast -o $@ estep_mem_resp.cc ../variant_fraction.cc ../../../util/loglib_header.cc -lgsl -lgslcblas -lm

estep_mem_resp_grid: estep_mem_resp_grid.cc ../variant_fraction.cc ../setting.hh ../../../util/loglib_header.hh ../../../util/enumtree_mutation_em_vf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -Ofast -o $@ estep_mem_resp_grid.cc ../variant_fraction.cc ../../../util/loglib_header.cc -lgsl -lgslcblas -lm

estep_mem_resp_grid_multinomial: estep_mem_resp_grid_multinomial.cc ../variant_fraction.cc ../setting.hh ../../../util/loglib_header.hh ../../../util/enumtree_mutation_em_multinomial.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -Ofast -o $@ estep_mem_resp_grid_multinomial.cc ../variant_fraction.cc ../../../util/loglib_header.cc -lgsl -lgslcblas -lm

em_with_normal: em_with_normal.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_mutation_em.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -Ofast -o $@ em_with_normal.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

accuracy: accuracy.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_wf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ accuracy.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

accuracy_bfgs2: accuracy_bfgs2.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_wf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ accuracy_bfgs2.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

calc_llik: calc_llik.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_wf_n_r_inherited.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -Ofast -o $@ calc_llik.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

calc_llik_u: calc_llik_u.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_wf_n_r_inherited.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -Ofast -o $@ calc_llik_u.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

growth: growth.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_wf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ growth.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

calc_params: calc_params.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_wf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ calc_params.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

calc_vf_dvf: calc_vf_dvf.cc ../variant_fraction.cc ../setting.hh ../../../util/loglib_header.hh ../../../util/enumtree_mutation_em_vf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -Ofast -o $@ calc_vf_dvf.cc ../variant_fraction.cc ../../../util/loglib_header.cc -lgsl -lgslcblas -lm

calc_vf_dvf_multinomial: calc_vf_dvf_multinomial.cc ../variant_fraction.cc ../setting.hh ../../../util/loglib_header.hh ../../../util/enumtree_mutation_em_multinomial.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -Ofast -o $@ calc_vf_dvf_multinomial.cc ../variant_fraction.cc ../../../util/loglib_header.cc -lgsl -lgslcblas -lm

calc_reverse_param: calc_reverse_param.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_mutation_em.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ calc_reverse_param.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

calc_gsl2param: calc_gsl2param.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_mutation_em.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ calc_gsl2param.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

growth_u_beta: growth_u_beta.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_wf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ growth_u_beta.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

%: %.cc setting.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ $< -lgsl -lgslcblas -lm

.PHONY: check-syntax

check-syntax:
	$(CXX) -Wall -Wextra -pedantic -fsyntax-only em.cc

