#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "stb_image.h"
#include "stb_image_write.h"
// https://github.com/nothings/stb
// biblioteca de read/write imagini fara optimizari

#include <bits/stdc++.h>

using namespace std;

const int BLOCK_SIZE = 8;

// taken from https://en.wikipedia.org/wiki/Quantization_(image_processing)
const int Q[8][8] = {
    {16, 11, 10, 16, 24, 40, 51, 61},
    {12, 12, 14, 19, 26, 58, 60, 55},
    {14, 13, 16, 24, 40, 57, 69, 56},
    {14, 17, 22, 29, 51, 87, 80, 62},
    {18, 22, 37, 56, 68, 109, 103, 77},
    {24, 35, 55, 64, 81, 104, 113, 92},
    {49, 64, 78, 87, 103, 121, 120, 101},
    {72, 92, 95, 98, 112, 100, 103, 99}};

const int zigzag_index[64][2] = {
    {0,0}, {0,1}, {1,0}, {2,0}, {1,1}, {0,2}, {0,3}, {1,2},
    {2,1}, {3,0}, {4,0}, {3,1}, {2,2}, {1,3}, {0,4}, {0,5},
    {1,4}, {2,3}, {3,2}, {4,1}, {5,0}, {6,0}, {5,1}, {4,2},
    {3,3}, {2,4}, {1,5}, {0,6}, {0,7}, {1,6}, {2,5}, {3,4},
    {4,3}, {5,2}, {6,1}, {7,0}, {7,1}, {6,2}, {5,3}, {4,4},
    {3,5}, {2,6}, {1,7}, {2,7}, {3,6}, {4,5}, {5,4}, {6,3},
    {7,2}, {7,3}, {6,4}, {5,5}, {4,6}, {3,7}, {4,7}, {5,6},
    {6,5}, {7,4}, {7,5}, {6,6}, {5,7}, {6,7}, {7,6}, {7,7}
};

void dct_block(const vector<vector<double>> &block, vector<vector<double>> &coeff)
{
    for (int u = 0; u < BLOCK_SIZE; ++u)
    {
        for (int v = 0; v < BLOCK_SIZE; ++v)
        {
            double sum = 0.0;
            for (int x = 0; x < BLOCK_SIZE; ++x)
            {
                for (int y = 0; y < BLOCK_SIZE; ++y)
                {
                    sum += block[x][y] * cos((2 * x + 1) * u * M_PI / 16.0) * cos((2 * y + 1) * v * M_PI / 16.0);
                }
            }
            double cu = (u == 0) ? 1.0 / sqrt(2.0) : 1.0;
            double cv = (v == 0) ? 1.0 / sqrt(2.0) : 1.0;
            coeff[u][v] = 0.25 * cu * cv * sum;
        }
    }
}

void idct_block(const vector<vector<double>> &coeff, vector<vector<double>> &block)
{
    for (int x = 0; x < BLOCK_SIZE; ++x)
    {
        for (int y = 0; y < BLOCK_SIZE; ++y)
        {
            double sum = 0.0;
            for (int u = 0; u < BLOCK_SIZE; ++u)
            {
                for (int v = 0; v < BLOCK_SIZE; ++v)
                {
                    double cu = (u == 0) ? 1.0 / sqrt(2.0) : 1.0;
                    double cv = (v == 0) ? 1.0 / sqrt(2.0) : 1.0;
                    sum += cu * cv * coeff[u][v] * cos((2 * x + 1) * u * M_PI / 16.0) * cos((2 * y + 1) * v * M_PI / 16.0);
                }
            }
            block[x][y] = 0.25 * sum;
        }
    }
}

void compress_channel(const vector<vector<double>> &channel, vector<vector<double>> &output, int width, int height, int K)
{
    output.resize(height, vector<double>(width, 0));

    // parcurgem imaginea in blocuri de 8x8
    for (int i = 0; i < height; i += BLOCK_SIZE)
    {
        for (int j = 0; j < width; j += BLOCK_SIZE)
        {
            vector<vector<double>> block(BLOCK_SIZE, vector<double>(BLOCK_SIZE, 0));
            vector<vector<double>> coeff(BLOCK_SIZE, vector<double>(BLOCK_SIZE, 0));

            // in caz ca imaginea nu e multiplu de 8
            int block_h = min(BLOCK_SIZE, height - i);
            int block_w = min(BLOCK_SIZE, width - j);

            // extragem blocul din canal
            for (int x = 0; x < block_h; ++x)
                for (int y = 0; y < block_w; ++y)
                    block[x][y] = channel[i + x][j + y];

            dct_block(block, coeff);

            // cuantizare (impartim fiecare coeficient la valoarea din matricea Q)
            for (int u = 0; u < BLOCK_SIZE; ++u)
                for (int v = 0; v < BLOCK_SIZE; ++v)
                    coeff[u][v] = round(coeff[u][v] / Q[u][v]);

            // convertim in vector 1D in ordine zigzag
            vector<double> zz(64);
            for (int idx = 0; idx < 64; ++idx)
            {
                int u = zigzag_index[idx][0];
                int v = zigzag_index[idx][1];
                zz[idx] = coeff[u][v];
            }

            // pastram doar primii K coeficienti
            for (int idx = K; idx < 64; ++idx)
                zz[idx] = 0;

            // reconvertim vectorul zigzag in matrice 8x8 si dequantizam
            for (int idx = 0; idx < 64; ++idx)
            {
                int u = zigzag_index[idx][0];
                int v = zigzag_index[idx][1];
                coeff[u][v] = zz[idx] * Q[u][v];
            }

            // reconstruim blocul
            idct_block(coeff, block);

            // punem blocul reconstruit in canalul de output
            for (int x = 0; x < block_h; ++x)
                for (int y = 0; y < block_w; ++y)
                    output[i + x][j + y] = block[x][y];
        }
    }
}

int main(int argc, char *argv[])
{
    string input_file;
    string output_file = "compressed_color.jpg";
    string output_original = "original_copy.jpg";
    int K;
    if (argc == 3)
    {
        if(atoi(argv[2]) < 1 || atoi(argv[2]) > 64) {
            cerr << "K must be between 1 and 64" << endl;
            return 1;
        } else {
            input_file = argv[1];
            K = atoi(argv[2]);
        }
    }
    else if (argc == 2)
    {
        input_file = argv[1];
        K = 5;
    }
    else
    {
        cerr << "Usage: " << argv[0] << " <input_image> [K (1-64, default 10)]" << endl;
        return 1;
    }

    

    int width, height, channels;
    unsigned char *img = stbi_load(input_file.c_str(), &width, &height, &channels, 0);
    if (!img)
    {
        cerr << "Image can't be read:" << input_file << endl;
        return 1;
    }

    cout << "Image successfully loaded " << width << "x" << height << ", Nr channels: " << channels << endl;

    // salvam copia originala ca sa putem sa facem comparatia fara optimizarile cu care vine fisierul jpeg original
    if (!stbi_write_jpg(output_original.c_str(), width, height, channels, img, 90))
    {
        cerr << "Error writing copy image" << endl;
        stbi_image_free(img);
        return 1;
    }

    vector<vector<double>> R(height, vector<double>(width, 0));
    vector<vector<double>> G(height, vector<double>(width, 0));
    vector<vector<double>> B(height, vector<double>(width, 0));

    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            int idx = (y * width + x) * channels;
            R[y][x] = img[idx];
            G[y][x] = (channels > 1) ? img[idx + 1] : R[y][x];
            B[y][x] = (channels > 2) ? img[idx + 2] : R[y][x];
        }
    }
    stbi_image_free(img);

    vector<vector<double>> R_comp, G_comp, B_comp;
    compress_channel(R, R_comp, width, height, K);
    compress_channel(G, G_comp, width, height, K);
    compress_channel(B, B_comp, width, height, K);

    vector<unsigned char> output(width * height * channels);
    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            int idx = (y * width + x) * channels;
            output[idx] = max(0.0, min(255.0, R_comp[y][x]));
            if (channels > 1)
                output[idx + 1] = max(0.0, min(255.0, G_comp[y][x]));
            if (channels > 2)
                output[idx + 2] = max(0.0, min(255.0, B_comp[y][x]));
        }
    }

    if (!stbi_write_jpg(output_file.c_str(), width, height, channels, output.data(), 90))
    {
        cerr << "Error writing compressed image" << endl;
        return 1;
    }

    cout << "Compressed image saved as: " << output_file << endl;
    cout << "Original image saved as: " << output_original << endl;
    return 0;
}
