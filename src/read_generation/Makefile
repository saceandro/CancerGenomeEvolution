CXX = g++

generate: generate.cc variant_fraction.cc setting.hh ../../util/loglib_header.hh ../../util/enumtree_wf_n_r.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ generate.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

generate_params: generate_params.cc variant_fraction.cc setting.hh ../../util/loglib_header.hh ../../util/enumtree_wf_n_r.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ generate_params.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

generate_reads: generate_reads.cc variant_fraction.cc setting.hh ../../util/loglib_header.hh ../../util/enumtree_wf_n_r.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ generate_reads.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

generate_reads_allinherited: generate_reads_allinherited.cc variant_fraction.cc setting.hh ../../util/loglib_header.hh ../../util/enumtree_wf_n_r.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -Ofast -o $@ generate_reads_allinherited.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

generate_reads_poisson: generate_reads_poisson.cc variant_fraction.cc setting.hh ../../util/loglib_header.hh ../../util/enumtree_wf_n_r.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -Ofast -o $@ generate_reads_poisson.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

variant_fraction: variant_fraction.cc variant_fraction.cc setting.hh ../../util/loglib_header.hh ../../util/enumtree_wf_n_r.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ variant_fraction.cc variant_fraction.cc ../../util/loglib_header.cc -lgsl -lgslcblas -lm

%: %.cc setting.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ $< -lgsl -lgslcblas -lm

.PHONY: check-syntax

check-syntax:
	$(CXX) -Wall -Wextra -pedantic -fsyntax-only em.cc

