#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "stb_image.h"
#include "stb_image_write.h"
// https://github.com/nothings/stb
// biblioteca de read/write imagini fara optimizari

#include <bits/stdc++.h>

using namespace std;

const int BLOCK_SIZE = 8;

// metode de imbunatatire: metoda zig-zag si pastrarea primilor K coeficienti din fiecare block

void dct_block(const vector<vector<double>>& block, vector<vector<double>>& coeff) {
    for (int u = 0; u < BLOCK_SIZE; ++u) {
        for (int v = 0; v < BLOCK_SIZE; ++v) {
            double sum = 0.0;
            for (int x = 0; x < BLOCK_SIZE; ++x) {
                for (int y = 0; y < BLOCK_SIZE; ++y) {
                    sum += block[x][y] *
                           cos((2*x+1)*u*M_PI/16.0) *
                           cos((2*y+1)*v*M_PI/16.0);
                }
            }
            double cu = (u == 0) ? 1.0/sqrt(2.0) : 1.0;
            double cv = (v == 0) ? 1.0/sqrt(2.0) : 1.0;
            coeff[u][v] = 0.25 * cu * cv * sum;
        }
    }
}

void idct_block(const vector<vector<double>>& coeff, vector<vector<double>>& block) {
    for (int x = 0; x < BLOCK_SIZE; ++x) {
        for (int y = 0; y < BLOCK_SIZE; ++y) {
            double sum = 0.0;
            for (int u = 0; u < BLOCK_SIZE; ++u) {
                for (int v = 0; v < BLOCK_SIZE; ++v) {
                    double cu = (u == 0) ? 1.0/sqrt(2.0) : 1.0;
                    double cv = (v == 0) ? 1.0/sqrt(2.0) : 1.0;
                    sum += cu * cv * coeff[u][v] *
                           cos((2*x+1)*u*M_PI/16.0) *
                           cos((2*y+1)*v*M_PI/16.0);
                }
            }
            block[x][y] = 0.25 * sum;
        }
    }
}

void compress_channel(const vector<vector<double>>& channel, vector<vector<double>>& output, int width, int height, double threshold) {
    output = channel;
    for (int i = 0; i < height; i += BLOCK_SIZE) {
        for (int j = 0; j < width; j += BLOCK_SIZE) {
            vector<vector<double>> block(BLOCK_SIZE, vector<double>(BLOCK_SIZE, 0));
            vector<vector<double>> coeff(BLOCK_SIZE, vector<double>(BLOCK_SIZE, 0));
            for (int x = 0; x < BLOCK_SIZE; ++x)
                for (int y = 0; y < BLOCK_SIZE; ++y)
                    if (i+x < height && j+y < width)
                        block[x][y] = channel[i+x][j+y];

            dct_block(block, coeff);

            for (int u = 0; u < BLOCK_SIZE; ++u)
                for (int v = 0; v < BLOCK_SIZE; ++v)
                    if (fabs(coeff[u][v]) < threshold)
                        coeff[u][v] = 0;

            idct_block(coeff, block);

            for (int x = 0; x < BLOCK_SIZE; ++x)
                for (int y = 0; y < BLOCK_SIZE; ++y)
                    if (i+x < height && j+y < width)
                        output[i+x][j+y] = block[x][y];
        }
    }
}

int main() {
    string input_file = "pisica512.jpg";
    string output_file = "compressed_color.jpg";
    string output_original = "original_copy.jpg";
    double threshold = 50.0;

    int width, height, channels;
    unsigned char* img = stbi_load(input_file.c_str(), &width, &height, &channels, 0);
    if (!img) {
        cerr << "Image can't be read:" << input_file << endl;
        return 1;
    }

    cout << "Image successfully loaded " << width << "x" << height << ", Nr channels: " << channels << endl;

    // salvam copia originala ca sa putem sa facem comparatia fara optimizarile cu care vine fisierul jpeg original
    if (!stbi_write_jpg(output_original.c_str(), width, height, channels, img, 90)) {
        cerr << "Error writing copy image" << endl;
        stbi_image_free(img);
        return 1;
    }

    vector<vector<double>> R(height, vector<double>(width, 0));
    vector<vector<double>> G(height, vector<double>(width, 0));
    vector<vector<double>> B(height, vector<double>(width, 0));

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            int idx = (y*width + x) * channels;
            R[y][x] = img[idx];
            G[y][x] = (channels > 1) ? img[idx+1] : R[y][x];
            B[y][x] = (channels > 2) ? img[idx+2] : R[y][x];
        }
    }
    stbi_image_free(img);

    vector<vector<double>> R_comp, G_comp, B_comp;
    compress_channel(R, R_comp, width, height, threshold);
    compress_channel(G, G_comp, width, height, threshold);
    compress_channel(B, B_comp, width, height, threshold);

    vector<unsigned char> output(width*height*channels);
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            int idx = (y*width + x)*channels;
            output[idx] = max(0.0, min(255.0, R_comp[y][x]));
            if (channels > 1) output[idx+1] = max(0.0, min(255.0, G_comp[y][x]));
            if (channels > 2) output[idx+2] = max(0.0, min(255.0, B_comp[y][x]));
        }
    }

    if (!stbi_write_jpg(output_file.c_str(), width, height, channels, output.data(), 90)) {
        cerr << "Error writing compressed image" << endl;
        return 1;
    }

    cout << "Compressed image saved as: " << output_file << endl;
    cout << "Original image saved as: " << output_original << endl;
    return 0;
}
