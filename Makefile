CFLAGS_      = $(CFLAGS) -O3 -std=gnu99 -Wall -Wno-unused
MCC_FLAGS_   = --ompss
MCC_FLAGS_I_ = --instrument
MCC_FLAGS_D_ = --debug -g -k
LDFLAGS_     = $(LDFLAGS)

## Check if coies are needed
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
	CFLAGS_  += -I$(MKL_INC_DIR) -DUSE_MKL
#	LDFLAGS_ += -L$(MKL_LIB_DIR) -lmkl_sequential -lmkl_core -lmkl_intel_lp64
	LDFLAGS_ += -L$(MKL_LIB_DIR) -lmkl_sequential -lmkl_core -lmkl_rt
else ifeq ($(OPENBLAS_SUPPORT_),YES)
	CFLAGS_  += -I$(OPENBLAS_INC_DIR) -DUSE_OPENBLAS
	LDFLAGS_ += -L$(OPENBLAS_LIB_DIR) -lopenblas
endif

PROGRAM_ = cholesky
PROGS_   = $(PROGRAM_)-p $(PROGRAM_)-i $(PROGRAM_)-d $(PROGRAM_)-seq

MCC  ?= mcc
MCC_  = $(CROSS_COMPILE)$(MCC)
GCC  ?= gcc
GCC_  = $(CROSS_COMPILE)$(GCC)

all: $(PROGS_)

$(PROGRAM_)-p: ./src/$(PROGRAM_).c
	$(MCC_) $(CFLAGS_) $(MCC_FLAGS_) $^ -o $@ $(LDFLAGS_)

$(PROGRAM_)-i:  ./src/$(PROGRAM_).c
	$(MCC_) $(CFLAGS_) $(MCC_FLAGS_) $(MCC_FLAGS_I_) $^ -o $@ $(LDFLAGS_)

$(PROGRAM_)-d:  ./src/$(PROGRAM_).c
	$(MCC_) $(CFLAGS_) $(MCC_FLAGS_) $(MCC_FLAGS_D_) $^ -o $@ $(LDFLAGS_)

$(PROGRAM_)-seq: ./src/$(PROGRAM_).c
	$(GCC_) $(CFLAGS_) $^ -o $@ $(LDFLAGS_)

clean:
	rm -f *.o $(PROGS_) $(MCC_)_$(PROGRAM_).c *_hls_automatic_mcxx.cpp *.config cholesky_vivado.log
	rm -rf cholesky_vivado
