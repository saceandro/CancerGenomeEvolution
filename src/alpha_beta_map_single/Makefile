CXX = g++

alpha_beta_map_single_compare_gradient: alpha_beta_map_single_compare_gradient.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_wf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ alpha_beta_map_single_compare_gradient.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

alpha_beta_map_single_grad: alpha_beta_map_single_grad.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_wf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ alpha_beta_map_single_grad.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

alpha_beta_map_single_calc_llik: alpha_beta_map_single_calc_llik.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_wf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ alpha_beta_map_single_calc_llik.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

alpha_beta_single_generate: alpha_beta_single_generate.cc variant_fraction.cc setting.hh ../../util/loglib_header.hh ../../util/loglib_header.cc ../../util/enumtree_wf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ alpha_beta_single_generate.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

alpha_beta_map_single_accuracy: alpha_beta_map_single_accuracy.cc variant_fraction.cc setting.hh ../../util/loglib_header.hh ../../util/loglib_header.cc ../../util/enumtree_wf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ alpha_beta_map_single_accuracy.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

%: %.cc setting.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ $< -lgsl -lgslcblas -lm

.PHONY: check-syntax

check-syntax:
	$(CXX) -Wall -Wextra -pedantic -fsyntax-only em.cc

