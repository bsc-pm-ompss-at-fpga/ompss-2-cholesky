.PHONY: clean info
all: help

PROGRAM_     = cholesky

#CLANG       ?= clang
#CLANG_FLAGS = -fompss-2 -DRUNTIME_MODE=\"perf\"
#MCC         ?= fpgacc
#MCC_         = $(CROSS_COMPILE)$(MCC)
#GCC_         = $(CROSS_COMPILE)gcc
#CFLAGS_      = $(CFLAGS) -O3 -std=gnu99 -Wall -Wno-unused -Wno-unknown-pragmas
#MCC_FLAGS_   = $(MCC_FLAGS) --ompss-2 --fpga -DRUNTIME_MODE=\"perf\" --variable=fpga_check_limits_memory_port:0 --fpga-link
#MCC_FLAGS_I_ = $(MCC_FLAGS_) --instrument -DRUNTIME_MODE=\"instr\"
#MCC_FLAGS_D_ = $(MCC_FLAGS_) --debug -g -k -DRUNTIME_MODE=\"debug\"
#LDFLAGS_     = $(LDFLAGS) -lm

# FPGA bitstream Variables
FPGA_HWRUNTIME         ?= pom
FPGA_CLOCK             ?= 200
FPGA_MEMORY_PORT_WIDTH ?= 128
INTERCONNECT_OPT       ?= performance
SYRK_NUM_ACCS          ?= 1
GEMM_NUM_ACCS          ?= 1
TRSM_NUM_ACCS          ?= 1
BLOCK_SIZE             ?= 32
POTRF_SMP              ?= 1
FPGA_GEMM_II           ?= 1
FPGA_OTHER_II          ?= 1

# Include the corresponding compiler makefile
--setup: FORCE
  ifeq ($(COMPILER),llvm)
    include llvm.mk
  else
    ifeq ($(COMPILER),mcxx)
      include mcxx.mk
    else
      $(info No valid COMPILER variable defined, using mcxx)
      include mcxx.mk
    endif
  endif
FORCE:

ifeq ($(POTRF_SMP),1)
	CFLAGS_ += -DPOTRF_SMP
endif

## Check if copies are needed
ifneq (,$(findstring USE_DMA_MEM,$(CFLAGS_)))
	# When using Kernel memory no copies are needed
	MCC_FLAGS_ += --no-copy-deps
else
	# Copy deps is the default
endif

## MKL Variables
MKL_DIR      ?= $(MKLROOT)
MKL_INC_DIR  ?= $(MKL_DIR)/include
MKL_LIB_DIR  ?= $(MKL_DIR)/lib
MKL_SUPPORT_ = $(if $(and $(wildcard $(MKL_INC_DIR)/mkl.h ), \
               $(wildcard $(MKL_LIB_DIR)/libmkl_sequential.so )),YES,NO)

## Open Blas Variables
OPENBLAS_DIR      ?= $(OPENBLAS_HOME)
OPENBLAS_INC_DIR  ?= $(OPENBLAS_DIR)/include
OPENBLAS_LIB_DIR  ?= $(OPENBLAS_DIR)/lib
OPENBLAS_SUPPORT_ = $(if $(and $(wildcard $(OPENBLAS_INC_DIR)/lapacke.h ), \
                    $(wildcard $(OPENBLAS_LIB_DIR)/libopenblas.so )),YES,NO)

ifeq ($(MKL_SUPPORT_),YES)
	COMPILER_FLAGS_  += -I$(MKL_INC_DIR) -DUSE_MKL
	LINKER_FLAGS_ += -L$(MKL_LIB_DIR) -lmkl_sequential -lmkl_core -lmkl_rt
else ifeq ($(OPENBLAS_SUPPORT_),YES)
	COMPILER_FLAGS_  += -I$(OPENBLAS_INC_DIR) -DUSE_OPENBLAS
	LINKER_FLAGS_ += -L$(OPENBLAS_LIB_DIR) -lopenblas -Wl,-rpath=$(OPENBLAS_LIB_DIR)
endif

COMPILER_FLAGS_ += -DFPGA_OTHER_LOOP_II=$(FPGA_OTHER_II) -DFPGA_GEMM_LOOP_II=$(FPGA_GEMM_II) -DBLOCK_SIZE=$(BLOCK_SIZE) -DFPGA_MEMORY_PORT_WIDTH=$(FPGA_MEMORY_PORT_WIDTH) -DSYRK_NUM_ACCS=$(SYRK_NUM_ACCS) -DGEMM_NUM_ACCS=$(GEMM_NUM_ACCS) -DTRSM_NUM_ACCS=$(TRSM_NUM_ACCS) -DBOARD=\"$(BOARD)\"

COMPILER_FLAGS_   += -DRUNTIME_MODE=\"perf\"
COMPILER_FLAGS_D_ += -DRUNTIME_MODE=\"debug\"
COMPILER_FLAGS_I_ += -DRUNTIME_MODE=\"instr\"

ifdef USE_URAM
	CFLAGS += -DUSE_URAM
endif

common-help:
	@echo 'Supported targets:       $(PROGRAM_)-p, $(PROGRAM_)-i, $(PROGRAM_)-d, $(PROGRAM_)-seq, design-p, design-i, design-d, bitstream-p, bitstream-i, bitstream-d, clean, help'
	@echo 'Environment variables:   CFLAGS, CROSS_COMPILE, LDFLAGS, MCC, MCC_FLAGS'
	@echo 'FPGA env. variables:     BOARD, FPGA_HWRUNTIME, FPGA_CLOCK, FPGA_MEMORY_PORT_WIDTH, SYRK_NUM_ACCS, GEMM_NUM_ACCS, TRSM_NUM_ACCS, BLOCK_SIZE, POTRF_SMP'
	@echo 'MKL env. variables:      MKLROOT, MKL_DIR, MKL_INC_DIR, MKL_LIB_DIR'
	@echo 'OpenBLAS env. variables: OPENBLAS_HOME, OPENBLAS_DIR, OPENBLAS_INC_DIR, OPENBLAS_LIB_DIR, OPENBLAS_IMPL'

$(PROGRAM_)-p: ./src/$(PROGRAM_).c
	$(COMPILER_) $(COMPILER_FLAGS_) $^ -o $@ $(LINKER_FLAGS_)

$(PROGRAM_)-i: ./src/$(PROGRAM_).c
	$(COMPILER_) $(COMPILER_FLAGS_) $(COMPILER_FLAGS_I_) $^ -o $@ $(LINKER_FLAGS_)

$(PROGRAM_)-d: ./src/$(PROGRAM_).c
	$(COMPILER_) $(COMPILER_FLAGS_) $(COMPILER_FLAGS_D_) $^ -o $@ $(LINKER_FLAGS_)

$(PROGRAM_)-seq: ./src/$(PROGRAM_).c
	$(COMPILER_) $(COMPILER_FLAGS_) $^ -o $@ $(LINKER_FLAGS_)

design-p: ./src/$(PROGRAM_).c
	$(eval TMPFILE := $(shell mktemp))
	$(COMPILER_) $(COMPILER_FLAGS_) \
		$(AIT_FLAGS_) $(AIT_FLAGS_DESIGN_) \
		$^ -o $(TMPFILE) $(LINKER_FLAGS_)
	rm $(TMPFILE)

design-i: ./src/$(PROGRAM_).c
	$(eval TMPFILE := $(shell mktemp))
	$(COMPILER_) $(COMPILER_FLAGS_I_) \
		$(AIT_FLAGS_) $(AIT_FLAGS_DESIGN_) \
		$^ -o $(TMPFILE) $(LINKER_FLAGS_)
	rm $(TMPFILE)

design-d: ./src/$(PROGRAM_).c
	$(eval TMPFILE := $(shell mktemp))
	$(COMPILER_) $(COMPILER_FLAGS_D_) \
		$(AIT_FLAGS_) $(AIT_FLAGS_DESIGN_) $(AIT_FLAGS_D_) \
		$^ -o $(TMPFILE) $(LINKER_FLAGS_)
	rm $(TMPFILE)

bitstream-p: ./src/$(PROGRAM_).c
	$(eval TMPFILE := $(shell mktemp))
	$(COMPILER_) $(COMPILER_FLAGS_) \
		$(AIT_FLAGS_) \
		$^ -o $(TMPFILE) $(LINKER_FLAGS_)
	rm $(TMPFILE)

bitstream-i: ./src/$(PROGRAM_).c
	$(eval TMPFILE := $(shell mktemp))
	$(COMPILER_) $(COMPILER_FLAGS_I_) \
		$(AIT_FLAGS_) \
		$^ -o $(TMPFILE) $(LINKER_FLAGS_)
	rm $(TMPFILE)

bitstream-d: ./src/$(PROGRAM_).c
	$(eval TMPFILE := $(shell mktemp))
	$(COMPILER_) $(COMPILER_FLAGS_D_) \
		$(AIT_FLAGS_) $(AIT_FLAGS_D_) \
		$^ -o $(TMPFILE) $(LINKER_FLAGS_)
	rm $(TMPFILE)

info:
	@echo "========== OPENBLAS =========="
	@echo "  SUPPORT enabled:  $(OPENBLAS_SUPPORT_)"
	@echo "  OPENBLAS_DIR:     $(OPENBLAS_DIR)"
	@echo "  OPENBLAS_INC_DIR: $(OPENBLAS_INC_DIR)"
	@echo "  OPENBLAS_LIB_DIR: $(OPENBLAS_LIB_DIR)"
	@echo "  Headers:          $(if $(wildcard $(OPENBLAS_INC_DIR)/lapacke.h ),YES,NO)"
	@echo "  Lib files (.so):  $(if $(wildcard $(OPENBLAS_LIB_DIR)/libopenblas.so ),YES,NO)"
	@echo "=============================="
	@echo "============= MKL ============"
	@echo "  SUPPORT enabled:  $(MKL_SUPPORT_)"
	@echo "  MKL_DIR:          $(MKL_DIR)"
	@echo "  MKL_INC_DIR:      $(MKL_INC_DIR)"
	@echo "  MKL_LIB_DIR:      $(MKL_LIB_DIR)"
	@echo "  Headers:          $(if $(wildcard $(MKL_INC_DIR)/mkl.h ),YES,NO)"
	@echo "  Lib files (.so):  $(if $(wildcard $(MKL_LIB_DIR)/libmkl_sequential.so ),YES,NO)"
	@echo "=============================="
