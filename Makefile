# RapidRefold Makefile
#
# Cross-platform build for Linux and macOS.
# Requires ViennaRNA submodule to be initialized and built first.
#
# Usage:
#   make setup          # Initialize submodule and build ViennaRNA (one-time)
#   make                # Build rapidrefold and tests
#   make test           # Run test suite
#   make clean          # Remove build artifacts
#   make distclean      # Also clean ViennaRNA build

CC       ?= gcc
CFLAGS   ?= -O2 -Wall -Wextra -Wno-unused-parameter
INCLUDES  = -I extern/ViennaRNA/src -I src

# ViennaRNA static library (built from submodule)
VRNA_LIB  = extern/ViennaRNA/src/ViennaRNA/.libs/libRNA_conv.a

# Platform-specific flags
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
  PLATFORM_LIBS = -lrt -lpthread
endif
ifeq ($(UNAME_S),Darwin)
  PLATFORM_LIBS = -lpthread
endif

# C++ runtime needed by ViennaRNA
ifeq ($(UNAME_S),Darwin)
  CXX_LIB = -lc++
else
  CXX_LIB = -lstdc++
endif

LDLIBS = $(VRNA_LIB) $(PLATFORM_LIBS) -lm $(CXX_LIB)

# Build directory
BUILD = build

# Source objects
OBJS = $(BUILD)/mfe_incremental.o

# Targets
.PHONY: all setup test clean distclean check_vrna

all: check_vrna $(BUILD)/rapidrefold $(BUILD)/fulltopk_experiment $(BUILD)/test_mfe

check_vrna:
	@if [ ! -f $(VRNA_LIB) ]; then \
	  echo ""; \
	  echo "Error: ViennaRNA library not found at $(VRNA_LIB)"; \
	  echo ""; \
	  echo "Run 'make setup' first to build ViennaRNA from the submodule."; \
	  echo ""; \
	  exit 1; \
	fi

# ---- Setup: initialize submodule and build ViennaRNA ----

setup:
	@echo "=== Checking build prerequisites ==="
	@missing=""; \
	command -v autoconf  >/dev/null 2>&1 || missing="$$missing autoconf"; \
	command -v automake  >/dev/null 2>&1 || missing="$$missing automake"; \
	(command -v libtoolize >/dev/null 2>&1 || command -v glibtoolize >/dev/null 2>&1) || missing="$$missing libtool"; \
	command -v $(CC)     >/dev/null 2>&1 || missing="$$missing cc"; \
	command -v xxd       >/dev/null 2>&1 || missing="$$missing xxd"; \
	if [ -n "$$missing" ]; then \
	  echo ""; \
	  echo "Error: missing required tools:$$missing"; \
	  echo ""; \
	  echo "Install them with:"; \
	  echo "  Ubuntu/Debian : sudo apt-get install autoconf automake libtool build-essential xxd"; \
	  echo "  Fedora/RHEL   : sudo dnf install autoconf automake libtool gcc gcc-c++ vim-common"; \
	  echo "  macOS         : brew install autoconf automake libtool"; \
	  echo ""; \
	  exit 1; \
	fi
	@echo "=== Initializing ViennaRNA submodule ==="
	git submodule update --init --recursive
	@echo ""
	@echo "=== Preparing ViennaRNA sources ==="
	@if [ -f extern/ViennaRNA/src/dlib-20.0.tar.bz2 ] && \
	    [ ! -d extern/ViennaRNA/src/dlib-20.0 ]; then \
	  echo "  Extracting dlib dependency..."; \
	  tar -xjf extern/ViennaRNA/src/dlib-20.0.tar.bz2 -C extern/ViennaRNA/src/; \
	fi
	@if ! autoconf --version 2>/dev/null | head -1 | awk '{split($$NF,v,"."); exit (v[1]>2||(v[1]==2&&v[2]>=71))?0:1}'; then \
	  echo "  Adding AC_C__GENERIC stub for autoconf < 2.71..."; \
	  printf 'dnl AC_C__GENERIC stub for autoconf < 2.71\nAC_DEFUN([AC_C__GENERIC], [])\n' \
	    > extern/ViennaRNA/m4/ac_c__generic.m4; \
	fi
	cd extern/ViennaRNA && autoreconf -fi 2>/dev/null || true
	@echo "  Patching compiler compatibility..."
	@sed -i.bak 's|double NAN = 1.0/0.0 - 1.0/0.0;|double NAN = __builtin_nan("");|' \
	  extern/ViennaRNA/src/cephes/const.c 2>/dev/null || true
	@sed -i.bak 's|double INFINITY = 1.0/0.0;|double INFINITY = __builtin_inf();|' \
	  extern/ViennaRNA/src/cephes/const.c 2>/dev/null || true
	@rm -f extern/ViennaRNA/src/cephes/const.c.bak
	@echo ""
	@echo "=== Configuring ViennaRNA ==="
	cd extern/ViennaRNA && \
	  ./configure --without-perl --without-python --without-swig \
	              --without-doc --without-tutorial --without-check \
	              --without-forester --without-cluster --without-kinwalker \
	              --without-json --without-cla \
	              --disable-shared --enable-static --disable-lto \
	              CFLAGS="-O2" 2>&1 | tail -5
	@echo ""
	@echo "=== Compiling ViennaRNA library (this may take a few minutes) ==="
	$(MAKE) -C extern/ViennaRNA/src/ViennaRNA -j$$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4) 2>&1 | tail -5
	@echo ""
	@if [ -f $(VRNA_LIB) ]; then \
	  echo "=== ViennaRNA built successfully ==="; \
	else \
	  echo "Error: ViennaRNA build failed. Run the following for full output:"; \
	  echo "  make -C extern/ViennaRNA/src/ViennaRNA"; \
	  exit 1; \
	fi

# ---- Build directory ----

$(BUILD):
	mkdir -p $(BUILD)

# ---- Object files ----

$(BUILD)/mfe_incremental.o: src/mfe_incremental.c src/mfe_incremental.h | $(BUILD)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# ---- Executables ----

$(BUILD)/rapidrefold: tools/rapidrefold.c $(OBJS) | $(BUILD)
	$(CC) $(CFLAGS) $(INCLUDES) -I src -o $@ $< $(OBJS) $(LDLIBS)

$(BUILD)/fulltopk_experiment: tools/incremental_fulltopk_experiment.c $(OBJS) | $(BUILD)
	$(CC) $(CFLAGS) $(INCLUDES) -I src -o $@ $< $(OBJS) $(LDLIBS)

$(BUILD)/test_mfe: tests/test_incremental_mfe.c $(OBJS) | $(BUILD)
	$(CC) $(CFLAGS) $(INCLUDES) -I src -o $@ $< $(OBJS) $(LDLIBS)

# ---- Test ----

test: $(BUILD)/test_mfe
	@echo "=== Running Tests ==="
	$(BUILD)/test_mfe

# ---- Clean ----

clean:
	rm -rf $(BUILD)

distclean: clean
	-$(MAKE) -C extern/ViennaRNA distclean 2>/dev/null || true
