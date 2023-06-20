# This mk file is intended to be used only by the common Makefile

help: common-help
	@echo 'Compiler env. variables:     CFLAGS, CROSS_COMPILE, LDFLAGS, MCC_FLAGS'

COMPILER_         = $(CROSS_COMPILE)fpgacc
COMPILER_FLAGS_   = $(CFLAGS) -O3 $(MCC_FLAGS) --ompss-2 --fpga
COMPILER_FLAGS_I_ = --instrument
COMPILER_FLAGS_D_ = --debug -g -k
LINKER_FLAGS_     = $(LDFLAGS)

AIT_FLAGS_        = --bitstream-generation --Wf,--name=$(PROGRAM_),--board=$(BOARD),-c=$(FPGA_CLOCK)
AIT_FLAGS_DESIGN_ = --Wf,--to_step=design
AIT_FLAGS_D_      = --Wf,--debug_intfs=both -k -i -v

#Picos config
AIT_FLAGS_ =--Wf,--max_deps_per_task=3,--max_args_per_task=11,--max_copies_per_task=11,--picos_tm_size=32,--picos_dm_size=102,--picos_vm_size=102

# Optional optimization FPGA variables
ifdef FPGA_MEMORY_PORT_WIDTH
	COMPILER_FLAGS_ += --variable=fpga_memory_port_width:$(FPGA_MEMORY_PORT_WIDTH)
endif
ifdef MEMORY_INTERLEAVING_STRIDE
	AIT_FLAGS_ += --Wf,--memory_interleaving_stride=$(MEMORY_INTERLEAVING_STRIDE)
endif
ifdef SIMPLIFY_INTERCONNECTION
	AIT_FLAGS_ += --Wf,--simplify_interconnection
endif
ifdef INTERCONNECT_PRIORITIES
	AIT_FLAGS_ += --Wf,--interconnect_priorities
endif
ifdef INTERCONNECT_OPT
	AIT_FLAGS_ += --Wf,--interconnect_opt=$(INTERCONNECT_OPT)
endif
ifdef INTERCONNECT_REGSLICE
	AIT_FLAGS_ += --Wf,--interconnect_regslice=$(INTERCONNECT_REGSLICE)
endif
ifdef FLOORPLANNING_CONSTR
	AIT_FLAGS_ += --Wf,--floorplanning_constr=$(FLOORPLANNING_CONSTR)
endif
ifdef SLR_SLICES
	AIT_FLAGS_ += --Wf,--slr_slices=$(SLR_SLICES)
endif
ifdef PLACEMENT_FILE
	AIT_FLAGS_ += --Wf,--placement_file=$(PLACEMENT_FILE)
endif
ifdef DISABLE_UTILIZATION_CHECK
	AIT_FLAGS_ += --Wf,--disable_utilization_check
endif

clean:
	rm -fv *.o $(PROGRAM_)-? $(COMPILER_)_$(PROGRAM_)*.c *hls_auto_mcxx.cpp ait_$(PROGRAM_)*.json
	rm -frv $(PROGRAM_)_ait
