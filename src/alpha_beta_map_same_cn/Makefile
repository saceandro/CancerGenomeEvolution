CXX = g++

alpha_beta_map_same_cn_compare_gradient: alpha_beta_map_same_cn_compare_gradient.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_wf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ alpha_beta_map_same_cn_compare_gradient.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

alpha_beta_map_same_cn_compare_gradient2: alpha_beta_map_same_cn_compare_gradient2.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_wf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ alpha_beta_map_same_cn_compare_gradient2.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

alpha_beta_map_same_cn_compare_gradient_check: alpha_beta_map_same_cn_compare_gradient_check.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_wf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ alpha_beta_map_same_cn_compare_gradient_check.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

alpha_beta_map_same_cn_compare_gradient_test_index_conversion: alpha_beta_map_same_cn_compare_gradient_test_index_conversion.cc variant_fraction.cc setting.hh ../../util/loglib.hh ../../util/enumtree_wf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ alpha_beta_map_same_cn_compare_gradient_test_index_conversion.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

alpha_beta_same_cn_generate: alpha_beta_same_cn_generate.cc variant_fraction.cc setting.hh ../../util/loglib_header.hh ../../util/loglib_header.cc ../../util/enumtree_wf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ alpha_beta_same_cn_generate.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

alpha_beta_map_same_cn_grad: alpha_beta_map_same_cn_grad.cc variant_fraction.cc setting.hh ../../util/loglib_header.hh ../../util/loglib_header.cc ../../util/enumtree_wf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ alpha_beta_map_same_cn_grad.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

alpha_beta_map_same_cn_grad_init: alpha_beta_map_same_cn_grad_init.cc variant_fraction.cc setting.hh ../../util/loglib_header.hh ../../util/loglib_header.cc ../../util/enumtree_wf.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ alpha_beta_map_same_cn_grad_init.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

%: %.cc setting.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ $< -lgsl -lgslcblas -lm

.PHONY: check-syntax

check-syntax:
	$(CXX) -Wall -Wextra -pedantic -fsyntax-only em.cc

