build:
	g++ -O0 -o dft_filter_serial dft_filter_serial.cpp

clean:
	rm -f dft_filter_serial
	rm -f serial_filtered_dft.jpg
	rm -f serial_original_copy.jpg
.PHONY: build clean

run: build
	./dft_filter_serial
