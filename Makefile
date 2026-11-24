CXX ?= clang++
CXXFLAGS ?= -O3 -std=c++17

all: upmr_gen upmr_qmc

upmr_gen: upmr_gen.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

upmr_qmc: upmr_qmc.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

# quick smoke test: spin-boson pair at beta=5 against ED (about a minute)
check: upmr_gen
	python3 gen_models.py sb --betas 5 --spec gen/sb.spec
	./upmr_gen --spec gen/sb.spec --beta 5 --nonlocal 1 --sweeps 4e5 --therm 2e4 --seed 1 | grep -E 'pairing|PHASE|RESULT'

clean:
	rm -f upmr_gen upmr_qmc

.PHONY: all check clean
