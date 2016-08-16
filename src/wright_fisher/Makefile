CXX = g++

%: %.cc grad3.hh setting.hh
	$(CXX) -I/usr/local/include -L/usr/local/lib -g -o $@ $< -lgsl -lgslcblas -lm
	# $(CXX) -I/usr/local/Cellar/boost/1.58.0/include/ -g -o $@ $< -lm

.PHONY: check-syntax

check-syntax:
	$(CXX) -Wall -Wextra -pedantic -fsyntax-only em.cc

