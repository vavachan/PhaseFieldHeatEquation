###########################################################################
## Makefile for OpenLB's automatic CSE workflow
TYPE = operator

# Directory containing the target files with dynamics type information
TARGET_DIR = targets/$(TYPE)

# List of targets to extract expressions
EXT_TARGET_LIST = target_$(TYPE).txt

# Define the directory containing the template for extractions
EXT_DIR = extractions/$(TYPE)

# Define the directory containing the .out files
EXP_DIR = expressions/$(TYPE)

# Define the directory for dumping the .cse.h files
CSE_DIR = expressions/$(TYPE)/include

# File for including the generated .cse.h files
INCLUDE_HEADER := $(CSE_DIR)/generated_cse.h

##########################################################################
.PHONY: target_list
target_list:
	bash unite_targets.sh $(TARGET_DIR) target_tmp.txt $(TYPE);
	bash format_targets_$(TYPE).sh target_tmp.txt $(EXT_TARGET_LIST) $(TYPE);
	rm target_tmp.txt;

##########################################################################
.PHONY: decompose_targets
decompose_targets: $(target_list)
	mkdir -p $(CSE_DIR)/; \
	bash decompose_targets.sh $(EXT_TARGET_LIST) $(EXT_DIR) $(INCLUDE_HEADER)

EXT_TARGETS := $(wildcard $(EXT_DIR)/*.txt)

EXT_OBJECTS := $(EXT_TARGETS:.txt=.out)

$(EXT_DIR)/%.out: $(EXT_DIR)/%.txt
	@mkdir -p $(EXT_DIR)/$*;
	@if [ ! -f "$(EXP_DIR)/$*.out" ]; then \
		content=$$(cat $<); \
		cp $(EXT_DIR)/extract.cpp.template $(EXT_DIR)/Makefile $(EXT_DIR)/$*/.; \
		python extraction.py "$$content" "$(EXT_DIR)/$*/extract.cpp.template" "$<" "$(EXT_DIR)/$*/" "$(EXP_DIR)/" "$(TYPE)"; \
	else \
		echo "Skip existing $(EXP_DIR)/$*.out."; \
	fi

.PHONY: check_targets
check_targets: $(EXT_OBJECTS)
	bash check_extraction.sh $(EXT_DIR) $(EXP_DIR)

##########################################################################
EXP_SOURCES := $(wildcard $(EXP_DIR)/*.out)

CSE_OBJECTS := $(EXP_SOURCES:.out=.cse.h)

$(EXP_DIR)/%.cse.h: $(EXP_DIR)/%.out
	@mkdir -p $(CSE_DIR);
	@if [ ! -f "$(CSE_DIR)/$*.cse.h" ]; then \
		python optimize.py "$<" "$(EXP_DIR)/$(TYPE).cse.template" "$(CSE_DIR)/$*.cse.h"; \
	else \
		echo "Skip existing $(CSE_DIR)/$*.cse.h."; \
	fi

.PHONY: check_optimized
check_optimized: $(CSE_OBJECTS)
	bash check_optimized.sh $(EXP_DIR) $(CSE_DIR)
##########################################################################

prepare_extraction: target_list decompose_targets

extraction: $(EXT_OBJECTS) check_targets

optimize_cse: $(CSE_OBJECTS) check_optimized

##########################################################################

.PHONY: clean-cse
clean-cse:
	rm -f $(CSE_DIR)/*.cse.h $(INCLUDE_HEADER)

.PHONY: clean-extraction
clean-extraction:
	rm -f $(EXP_DIR)/*.out $(EXT_DIR)/*.txt
	@for dir in $(EXT_DIR)/*; do\
		if [ -d "$$dir" ]; then \
			rm -rf "$$dir"/; \
		fi \
	done

.PHONY: clean-targets
clean-targets:
	rm -f $(EXT_TARGET_LIST)

.PHONY: clean
clean: clean-extraction

.PHONY: clean-all
clean-all: clean-cse clean-extraction clean-targets

