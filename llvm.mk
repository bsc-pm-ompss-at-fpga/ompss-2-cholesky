# This mk file is intended to be used only by the common Makefile

help: common-help
	@echo 'Environment variables:    CFLAGS, CROSS_COMPILE, LDFLAGS'

COMPILER_         = clang
COMPILER_FLAGS_   = $(CFLAGS) -fompss-2 -fompss-fpga-wrapper-code
COMPILER_FLAGS_D_ = $(COMPILER_FLAGS_) -g -fompss-fpga-hls-tasks-dir $(PWD)
LINKER_FLAGS_     = $(LDFLAGS)

AIT_FLAGS__        = --name=$(PROGRAM_) --board=$(BOARD) -c=$(FPGA_CLOCK) -k
AIT_FLAGS_DESIGN__ = --to_step=design
AIT_FLAGS_D__      = --debug_intfs=both

# Optional optimization FPGA variables
ifdef FPGA_MEMORY_PORT_WIDTH
	COMPILER_FLAGS_ += -fompss-fpga-memory-port-width $(FPGA_MEMORY_PORT_WIDTH)
endif
ifdef MEMORY_INTERLEAVING_STRIDE
	AIT_FLAGS__ += --memory_interleaving_stride=$(MEMORY_INTERLEAVING_STRIDE)
endif
ifdef SIMPLIFY_INTERCONNECTION
	AIT_FLAGS__ += --simplify_interconnection
endif
ifdef INTERCONNECT_OPT
	AIT_FLAGS__ += --interconnect_opt=$(INTERCONNECT_OPT)
endif
ifdef INTERCONNECT_REGSLICE
	AIT_FLAGS__ += --interconnect_regslice=$(INTERCONNECT_REGSLICE)
endif
ifdef FLOORPLANNING_CONSTR
	AIT_FLAGS__ += --floorplanning_constr=$(FLOORPLANNING_CONSTR)
endif
ifdef SLR_SLICES
	AIT_FLAGS__ += --slr_slices=$(SLR_SLICES)
endif
ifdef PLACEMENT_FILE
	AIT_FLAGS__ += --placement_file=$(PLACEMENT_FILE)
endif

AIT_FLAGS_        = -fompss-fpga-ait-flags "$(AIT_FLAGS__)"
AIT_FLAGS_DESIGN_ = -fompss-fpga-ait-flags "$(AIT_FLAGS_DESIGN__)"
AIT_FLAGS_D_      = -fompss-fpga-ait-flags "$(AIT_FLAGS_D__)"

clean:
	rm -fv *.o $(PROGRAM_)-? $(PROGRAM_)_hls_automatic_clang.cpp ait_extracted.json
	rm -frv $(PROGRAM_)_ait
