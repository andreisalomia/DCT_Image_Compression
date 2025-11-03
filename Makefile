build:
	g++ -o dct_compress dct_image_compress.cpp

clean:
	rm -f dct_compress
.PHONY: build clean

run: build
	./dct_compress