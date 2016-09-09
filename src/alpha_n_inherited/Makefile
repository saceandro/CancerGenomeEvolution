CXX = g++

generate: generate.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_wf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ generate.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

generate_params: generate_params.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_wf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ generate_params.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

generate_reads: generate_reads.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_wf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ generate_reads.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

compare: compare.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_wf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ compare.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

compare_simple: compare_simple.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_wf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ compare_simple.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

grad: grad.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_wf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ grad.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

grad_foraccuracy: grad_foraccuracy.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_wf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ grad_foraccuracy.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

accuracy: accuracy.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_wf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ accuracy.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

accuracy_bfgs2: accuracy_bfgs2.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_wf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ accuracy_bfgs2.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

llik: llik.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_wf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ llik.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

growth: growth.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_wf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ growth.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

calc_params: calc_params.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_wf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ calc_params.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

growth_u_beta: growth_u_beta.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_wf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ growth_u_beta.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

%: %.cc setting.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ $< -lgsl -lgslcblas -lm

.PHONY: check-syntax

check-syntax:
	$(CXX) -Wall -Wextra -pedantic -fsyntax-only em.cc

