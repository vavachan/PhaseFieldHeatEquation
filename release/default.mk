# Include config.mk environment (optional)
-include $(OLB_ROOT)/config.mk

# Select mixed compilation mode if separate CUDA compiler is given
ifdef CUDA_CXX
include $(OLB_ROOT)/default.mixed.mk
# otherwise use single CXX for all of OpenLB
else
include $(OLB_ROOT)/default.single.mk
endif
