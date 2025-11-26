#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "../stb_image.h"
#include "../stb_image_write.h"
// https://github.com/nothings/stb
// biblioteca de read/write imagini fara optimizari

#include <bits/stdc++.h>
#include <complex>

#define THRESHOLD_FACTOR 0.001

using namespace std;

const double PI = M_PI;

vector<complex<double>> dft_1d(const vector<complex<double>> &signal)
{
    int N = signal.size();
    vector<complex<double>> result(N, 0);
    for (int k = 0; k < N; ++k)
    {
        for (int n = 0; n < N; ++n)
        {
            double angle = 2.0 * PI * k * n / N;
            complex<double> expo(cos(angle), -sin(angle));
            result[k] += signal[n] * expo;
        }
    }
    return result;
}

vector<complex<double>> idft_1d(const vector<complex<double>> &signal)
{
    int N = signal.size();
    vector<complex<double>> result(N, 0);
    for (int k = 0; k < N; ++k)
    {
        for (int n = 0; n < N; ++n)
        {
            double angle = 2.0 * PI * k * n / N;
            complex<double> expo(cos(angle), sin(angle));
            result[k] += signal[n] * expo;
        }
        result[k] /= N;
    }
    return result;
}

vector<vector<complex<double>>> dft_2d(const vector<vector<complex<double>>> &image)
{
    int rows = image.size();
    int cols = image[0].size();

    vector<vector<complex<double>>> temp(rows, vector<complex<double>>(cols));
    for (int i = 0; i < rows; ++i)
    {
        temp[i] = dft_1d(image[i]);
    }

    vector<vector<complex<double>>> result(rows, vector<complex<double>>(cols));
    for (int j = 0; j < cols; ++j)
    {
        vector<complex<double>> column(rows);
        for (int i = 0; i < rows; ++i)
        {
            column[i] = temp[i][j];
        }
        vector<complex<double>> transformed_col = dft_1d(column);
        for (int i = 0; i < rows; ++i)
        {
            result[i][j] = transformed_col[i];
        }
    }

    return result;
}

vector<vector<complex<double>>> idft_2d(const vector<vector<complex<double>>> &freq)
{
    int rows = freq.size();
    int cols = freq[0].size();

    vector<vector<complex<double>>> temp(rows, vector<complex<double>>(cols));
    for (int i = 0; i < rows; ++i)
    {
        temp[i] = idft_1d(freq[i]);
    }

    vector<vector<complex<double>>> result(rows, vector<complex<double>>(cols));
    for (int j = 0; j < cols; ++j)
    {
        vector<complex<double>> column(rows);
        for (int i = 0; i < rows; ++i)
        {
            column[i] = temp[i][j];
        }
        vector<complex<double>> transformed_col = idft_1d(column);
        for (int i = 0; i < rows; ++i)
        {
            result[i][j] = transformed_col[i];
        }
    }

    return result;
}

void compress_channel_fft(const vector<vector<double>> &channel, vector<vector<double>> &output, int width, int height)
{
    output.resize(height, vector<double>(width, 0));

    vector<vector<complex<double>>> complex_channel(height, vector<complex<double>>(width));
    for (int i = 0; i < height; ++i)
    {
        for (int j = 0; j < width; ++j)
        {
            complex_channel[i][j] = complex<double>(channel[i][j], 0);
        }
    }

    vector<vector<complex<double>>> freq_domain = dft_2d(complex_channel);

    double max_magnitude = 0.0;
    for (int i = 0; i < height; ++i)
    {
        for (int j = 0; j < width; ++j)
        {
            double mag = abs(freq_domain[i][j]);
            if (mag > max_magnitude)
            {
                max_magnitude = mag;
            }
        }
    }

    double threshold = max_magnitude * THRESHOLD_FACTOR;

    cout << "  Max magnitude: " << max_magnitude << ", Threshold: " << threshold << endl;

    int coeffs_kept = 0;
    for (int i = 0; i < height; ++i)
    {
        for (int j = 0; j < width; ++j)
        {
            if (abs(freq_domain[i][j]) < threshold)
            {
                freq_domain[i][j] = complex<double>(0, 0);
            }
            else
            {
                coeffs_kept++;
            }
        }
    }

    cout << "  Coefficients kept: " << coeffs_kept << "/" << (width * height) << " (" << (100.0 * coeffs_kept / (width * height)) << "%)" << endl;

    vector<vector<complex<double>>> reconstructed = idft_2d(freq_domain);

    for (int i = 0; i < height; ++i)
    {
        for (int j = 0; j < width; ++j)
        {
            output[i][j] = reconstructed[i][j].real();
        }
    }
}

int main(int argc, char *argv[])
{
    string input_file;
    string output_file = "serial_filtered_dft.jpg";
    string output_original = "serial_original_copy.jpg";

    if (argc == 2)
    {
        input_file = argv[1];
    }
    else
    {
        cerr << "Usage: " << argv[0] << " <input_image>" << endl;
        return 1;
    }

    int width, height, channels;
    unsigned char *img = stbi_load(input_file.c_str(), &width, &height, &channels, 0);
    if (!img)
    {
        cerr << "Image can't be read: " << input_file << endl;
        return 1;
    }

    cout << "Image successfully loaded " << width << "x" << height << ", Nr channels: " << channels << endl;

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

    cout << "Processing Red channel..." << endl;
    vector<vector<double>> R_comp, G_comp, B_comp;
    compress_channel_fft(R, R_comp, width, height);

    cout << "Processing Green channel..." << endl;
    compress_channel_fft(G, G_comp, width, height);

    cout << "Processing Blue channel..." << endl;
    compress_channel_fft(B, B_comp, width, height);

    vector<unsigned char> output(width * height * channels);
    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            int idx = (y * width + x) * channels;
            output[idx] = static_cast<unsigned char>(max(0.0, min(255.0, R_comp[y][x])));
            if (channels > 1)
                output[idx + 1] = static_cast<unsigned char>(max(0.0, min(255.0, G_comp[y][x])));
            if (channels > 2)
                output[idx + 2] = static_cast<unsigned char>(max(0.0, min(255.0, B_comp[y][x])));
        }
    }

    if (!stbi_write_jpg(output_file.c_str(), width, height, channels, output.data(), 90))
    {
        cerr << "Error writing filtered image" << endl;
        return 1;
    }

    cout << "Filtered image saved as: " << output_file << endl;
    cout << "Original image saved as: " << output_original << endl;
    return 0;
}
