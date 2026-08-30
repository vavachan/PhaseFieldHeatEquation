include $(OLB_ROOT)/rules.mk
all: dependencies core $(EXAMPLE)

INCLUDE_DIRS := $(OLB_ROOT)/src

ifneq ($(USE_EMBEDDED_DEPENDENCIES), OFF)
INCLUDE_DIRS += \
	$(OLB_ROOT)/external/zlib \
	$(OLB_ROOT)/external/tinyxml2

LDFLAGS := -L$(OLB_ROOT)/external/lib $(LDFLAGS)

EXTENSION := $(if $(filter EMSCRIPTEN,$(FEATURES)),.js,)

dependencies:
	$(MAKE) -C $(OLB_ROOT)/external

clean-dependencies:
	$(MAKE) -C $(OLB_ROOT)/external clean
else
.PHONY: dependencies clean-dependencies
endif

LDFLAGS += -L$(OLB_ROOT)/build/lib

core:
	$(MAKE) -C $(OLB_ROOT) core

CPP_FILES := $(wildcard *.cpp)
OBJ_FILES := $(CPP_FILES:.cpp=.o)
DEP_FILES := $(CPP_FILES:.cpp=.d)
WASM_FILES := $(CPP_FILES:.cpp=.wasm)
JS_FILES := $(CPP_FILES:.cpp=.js)

INCLUDE_FLAGS := -I$(subst $(EMPTY) $(EMPTY), -I,$(INCLUDE_DIRS))

%.d: $(CPP_FILES)
	@$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) $< -MM -MT $(@:.d=.o) >$@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) -c -o $@ $<

$(EXAMPLE): $(OBJ_FILES)
	$(CXX) $(OBJ_FILES) -o $@$(EXTENSION) -lolbcore $(LDFLAGS)

.PHONY: onlysample
onlysample: $(EXAMPLE)

.PHONY: run
run: $(EXAMPLE)
	./$(EXAMPLE)

.PHONY: clean-tmp
clean-tmp:
	rm -f tmp/*.* tmp/vtkData/*.* tmp/vtkData/data/*.* tmp/imageData/*.* tmp/imageData/data/*.* tmp/gnuplotData/*.* tmp/gnuplotData/data/*.*

.PHONY: clean-core
clean-core:
	make -C $(OLB_ROOT) clean-core

.PHONY: clean
clean: clean-tmp
	rm -f $(OBJ_FILES) $(DEP_FILES) $(WASM_FILES) $(JS_FILES) $(EXAMPLE)

-include $(DEP_FILES)
