IN_FILES   = boo.in foo.in
OUT_FILES = $(basename $(IN_FILES))
types = X Y
list = $(foreach a, $(types), $(foreach b, $(OUT_FILES),$(a)_$(b).out))
exe:
	@echo "$(list)"	
