CXX ?= g++

# ===== Common flags =====
CXXFLAGS = -std=c++14 -Wall -O2
LDFLAGS  =
LDLIBS   = -lpthread

# OpenMP flags (for IC and any target that uses omp pragmas)
OMPFLAGS = -fopenmp

# Sources (your project includes .cpp from main.cpp, so listing is for dependency tracking)
SRCS = src/main.cpp src/mygraph.cpp src/algs.cpp src/obj.cpp

# Default
all: run_dir maxcut revmax_nm revmax maxcov ic preproc preproc_ic

run_dir:
	@mkdir -p run

# ===== Main targets =====
maxcov: run_dir $(SRCS)
	$(CXX) src/main.cpp -o run/maxcov $(CXXFLAGS) $(LDFLAGS) -DMAXCOV $(LDLIBS)

maxcut: run_dir $(SRCS)
	$(CXX) src/main.cpp -o run/maxcut $(CXXFLAGS) $(LDFLAGS) -DMAXCUT $(LDLIBS)

revmax: run_dir $(SRCS)
	$(CXX) src/main.cpp -o run/revmax $(CXXFLAGS) $(LDFLAGS) -DREVMAX_MON $(LDLIBS)

revmax_nm: run_dir $(SRCS)
	$(CXX) src/main.cpp -o run/revmax_nm $(CXXFLAGS) $(LDFLAGS) -DREVMAX_NM $(LDLIBS)

# ===== IC (OpenMP) =====
ic: run_dir $(SRCS)
	$(CXX) src/main.cpp -o run/ic $(CXXFLAGS) $(OMPFLAGS) $(LDFLAGS) -DIC $(LDLIBS)

# ===== Debug targets =====
debug: run_dir $(SRCS)
	$(CXX) src/main.cpp -o run/maxcov_debug -std=c++14 -Wall -Og -g -DMAXCOV $(LDLIBS)

maxcut_debug: run_dir $(SRCS)
	$(CXX) src/main.cpp -o run/maxcut_debug -std=c++14 -Wall -Og -g -DMAXCUT $(LDLIBS)

revmax_debug: run_dir $(SRCS)
	$(CXX) src/main.cpp -o run/revmax_debug -std=c++14 -Wall -Og -g -DREVMAX_MON $(LDLIBS)

ic_debug: run_dir $(SRCS)
	$(CXX) src/main.cpp -o run/ic_debug -std=c++14 -Wall -Og -g -DIC $(OMPFLAGS) $(LDLIBS)

# ===== Tools =====
preproc: run_dir src/preprocess.cpp src/mygraph.cpp
	$(CXX) src/preprocess.cpp -o run/preproc $(CXXFLAGS) $(LDFLAGS) $(LDLIBS)

# ---- NEW: preprocess IC tool ----
# Place your file at src/preprocess_ic.cpp
# If you put it under src/data/preprocess_ic.cpp, change path accordingly.
preproc_ic: run_dir src/preprocess_ic.cpp src/mygraph.cpp
	$(CXX) src/preprocess_ic.cpp -o run/preproc_ic $(CXXFLAGS) $(LDFLAGS) $(LDLIBS)

er: run_dir src/gen_er.cpp src/mygraph.cpp
	$(CXX) -std=c++14 src/gen_er.cpp -o run/er -O2 -Wall

ba: run_dir src/gen_ba.cpp src/mygraph.cpp
	$(CXX) -std=c++14 src/gen_ba.cpp -o run/ba -O2 -Wall

clean:
	rm -f run/preproc run/preproc_ic run/maxcov run/revmax run/revmax_nm run/maxcut run/ic
	rm -f run/maxcov_debug run/maxcut_debug run/revmax_debug run/ic_debug

single-pass: maxcut preproc
	cd exp/$@; make DATA=er BIN=maxcut; make DATA=er BIN=revmax; \
	make DATA=ba BIN=maxcut; make DATA=ba BIN=revmax; \
	make DATA=fb BIN=maxcut; make DATA=fb BIN=revmax; \
	make DATA=slashdot BIN=maxcut; make DATA=slashdot BIN=revmax; \
	make DATA=pokec BIN=maxcut; make DATA=pokec BIN=revmax
