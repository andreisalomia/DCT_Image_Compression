build:
	g++ -o dft_compression dft_compression_serial.cpp

clean:
	rm -f dft_compression
.PHONY: build clean

run: build
	./dft_compression