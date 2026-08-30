all: lib zlib tinyxml2
.PHONY: zlib tinyxml2

lib:
	mkdir -p lib

zlib:
	make -C zlib
	cp zlib/build/libz.a lib/

clean_zlib:
	make -C zlib clean

tinyxml2:
	make -C tinyxml2
	cp tinyxml2/build/libtinyxml2.a lib/

clean_tinyxml2:
	make -C tinyxml2 clean

clean: clean_zlib clean_tinyxml2
	rm -f lib/libz.a lib/libtinyxml2.a