// condec.cpp

#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
#include <cstring>
#include <sstream>
#include <stdexcept>
#include <cmath>

using namespace std;

const double PI = 3.14159265358979323846;

// helper functions
struct WAVHeader {
    char riff[4];                // RIFF Header
    uint32_t chunkSize;           // RIFF Chunk Size
    char wave[4];                 // WAVE Header
    char fmt[4];                  // fmt subchunk
    uint32_t subchunk1Size;       // Size of the fmt chunk
    uint16_t audioFormat;         // Audio format (1=PCM, 3=IEEE float)
    uint16_t numChannels;         // Number of channels
    uint32_t sampleRate;          // Sampling rate
    uint32_t byteRate;            // (sampleRate * numChannels * bitsPerSample) / 8
    uint16_t blockAlign;          // numChannels * bitsPerSample / 8
    uint16_t bitsPerSample;       // Bits per sample
    char data[4];                 // "data" string
    uint32_t dataSize;            // Size of the data section
};
std::vector<std::vector<double>> read_wav(const char* filename, WAVHeader& header) {
    std::stringstream err;
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        err << "cannot open WAV file " << filename;
        throw std::runtime_error(err.str());
    }

    file.read(reinterpret_cast<char*>(&header), sizeof(WAVHeader));

    if (header.audioFormat != 1 && header.audioFormat != 3) { // Only PCM (1) and IEEE float (3) supported
        err << "unsupported WAV format (not PCM or IEEE float)";
        throw std::runtime_error(err.str());
    }

    int numSamples = header.dataSize / (header.bitsPerSample / 8);
    numSamples /= header.numChannels;
    std::vector<std::vector<double>> channels(header.numChannels, std::vector<double>(numSamples));

    for (int i = 0; i < numSamples; ++i) {
        for (int ch = 0; ch < header.numChannels; ++ch) {
            if (header.bitsPerSample == 16) {
                int16_t sample;
                file.read(reinterpret_cast<char*>(&sample), sizeof(int16_t));
                channels[ch][i] = sample / 32768.0;  // Normalize to [-1.0, 1.0]
            } else if (header.bitsPerSample == 32 && header.audioFormat == 3) {
                float sample;
                file.read(reinterpret_cast<char*>(&sample), sizeof(float));
                channels[ch][i] = sample;  // No need to normalize, as it's already in [-1.0, 1.0]
            } else {
                err << "unsupported bits per sample: " << header.bitsPerSample;
                throw std::runtime_error(err.str());
            }
        }
    }

    file.close();
    return channels;
}
void write_wav(const char* filename, std::vector<std::vector<double>>& channels, WAVHeader& header) {
    std::stringstream err;
    std::ofstream file(filename, std::ios::binary);

    if (!file) {
        err << "cannot write WAV file " << filename;
        throw std::runtime_error(err.str());
    }

    header.dataSize = 0;
    for (const auto& channel : channels) {
        header.dataSize += channel.size() * (header.bitsPerSample / 8);
    }
    // header.dataSize /= header.numChannels;
    file.write(reinterpret_cast<char*>(&header), sizeof(WAVHeader));
    int numSamples = channels[0].size();
    for (int i = 0; i < numSamples; ++i) {
        for (size_t ch = 0; ch < channels.size(); ++ch) {
            if (header.bitsPerSample == 16) {
                int16_t intSample = static_cast<int16_t>(channels[ch][i] * 32767);
                file.write(reinterpret_cast<char*>(&intSample), sizeof(int16_t));
            } else if (header.bitsPerSample == 32 && header.audioFormat == 3) {
                float floatSample = static_cast<float>(channels[ch][i]);
                file.write(reinterpret_cast<char*>(&floatSample), sizeof(float));
            } else {
                err << "unsupported bits per sample or audio format";
                throw std::runtime_error(err.str());
            }
        }
    }

    file.close();
}
std::string remove_extension(const std::string& filename) {
    size_t last_dot = filename.find_last_of(".");
    if (last_dot != std::string::npos) {
        return filename.substr(0, last_dot);
    }
    return filename; // No extension found, return the original filename
}

// math
int next_power_of_two (int n) {
    return pow (2, ceil (log2 (n)));
}
void fft (std::complex<double>* signal, int N) {
    if (N <= 1) return;

    std::vector<std::complex<double>> even (N / 2);
    std::vector<std::complex<double>> odd (N / 2);
    for (int i = 0; i < N / 2; ++i) {
        even[i] = signal[i * 2];
        odd[i] = signal[i * 2 + 1];
    }

    fft (even.data (), N / 2);
    fft (odd.data (), N / 2);

    for (int k = 0; k < N / 2; ++k) {
        std::complex<double> t = std::polar (1.0, -2 * PI * k / N) * odd[k];
        signal[k] = even[k] + t;
        signal[k + N / 2] = even[k] - t;
    }
}
void ifft (std::complex<double>* signal, int N) {
    for (int i = 0; i < N; ++i) {
        signal[i] = std::conj (signal[i]);
    }

    fft (signal, N);

    for (int i = 0; i < N; ++i) {
        signal[i] = std::conj (signal[i]) / static_cast<double> (N);
    }
}

double max_val_cplx (std::complex<double>* arr, int N, int& max_pos) {
    double max_magnitude = 0.0;
    max_pos = 0;  // Default position to 0

    for (int i = 0; i < N; ++i) {
        double magnitude = std::abs(arr[i]);
        if (magnitude > max_magnitude) {
            max_magnitude = magnitude;
            max_pos = i;
        }
    }

    return max_magnitude;
}
std::complex<double>* process (std::complex<double>* in_spectrum, std::complex<double>* deconv_signal, 
    std::complex<double>* result_spectrum, int N, double threshold, int op) {

    int mpos = 0;
    double mm = max_val_cplx (deconv_signal, N, mpos);

    for (int i = 0; i < N; ++i) {
        double denom_mag = std::abs (deconv_signal[i]);
        if (denom_mag > threshold * mm) {
            result_spectrum[i] = op >= 0 ? in_spectrum[i] * deconv_signal[i] : in_spectrum[i] / deconv_signal[i];
        } else {
            result_spectrum[i] = 0.0;
        }
    }

    return result_spectrum;
}

// ------------------------------
int main (int argc, char** argv) {
    cout << "[condec, ver. 0.1]" << endl << endl;
    cout << "(c) 2024 Carmine-Emanuele Cella" << endl << endl;
    try {
        if (argc != 8) {
            stringstream err;
            err << "syntax is '" << argv[0] << " <x_wav> <y_wav> <output_wav> op threshold scale mix'" << endl << endl
                << "where: " << endl
                << "op > 0 for convolution; op < 0 for deconvolution" << endl
                << "threshold is used for filtering in the frequency domain" << endl
                << "scale changes the amplification for <output_wav>" << endl
                << "mix is used to sum <x_wav> in conv or remove <y_wav> in deconv (residual)" << endl;

            throw runtime_error (err.str ());
        }

        const char* x_file = argv[1];
        const char* y_file = argv[2];
        const char* output_file = argv[3];
        const int op = atol (argv[4]);
        double threshold = atof (argv[5]);
        double scale = atof (argv[6]);
        double mix = atof (argv[7]);

        cout << "operation = " << (op > 0 ? "conv" : "deconv") << endl;
        cout << "threshold = " << threshold << endl;
        cout << "scale     = " << scale << endl;
        cout << "mix       = " << mix << endl;

        WAVHeader x_header, y_header;
        vector<vector<double>> x_data = read_wav (x_file, x_header);
        vector<vector<double>> y_data = read_wav (y_file, y_header);

        int max_ch = std::max (x_data.size (), y_data.size ());          
        int max_size = std::max (x_data[0].size(), y_data[0].size());
        int N = next_power_of_two (max_size);
        complex<double>* x_FFT = new complex<double>[N];
        complex<double>* y_FFT = new complex<double>[N];
        complex<double>* result_FFT = new complex<double>[N];
        
        int out_sz = N; //op >= 0 ? x_data[0].size () + y_data[0].size () : N;
        vector<vector<double>> output_data(max_ch, vector<double> (out_sz));
        vector<int> maxes (max_ch);

        cout << endl << "processing..."; cout.flush ();
        for (int i = 0; i < max_ch; ++i) {
            int ich = i >= (int) x_data.size () ? x_data.size () - 1 : i;
            int dch = i >= (int)  y_data.size () ? y_data.size () - 1 : i;
            x_data[ich].resize (N, 0.0);
            y_data[dch].resize (N, 0.0);

            for (int i = 0; i < N; ++i) {
                x_FFT[i] = x_data[ich][i];
                y_FFT[i] = y_data[dch][i];
            }

            fft (x_FFT, N);
            fft (y_FFT, N);
            process (x_FFT, y_FFT, result_FFT, N, threshold, op);
            ifft (result_FFT, N);

            max_val_cplx (result_FFT, N,  maxes[i]);
            for (int s = 0; s <  out_sz; ++s) {
                output_data[i][s] = result_FFT[s].real() * scale;
                if (s < (int) x_data[0].size () && op >= 0) output_data[i][s] += x_data[ich][s] * mix;
            }
        }
        cout << "done" << endl << endl;

        WAVHeader output_header = x_header;
        output_header.numChannels = output_data.size ();
        write_wav (output_file, output_data, output_header);
        cout << "output saved to  : " << output_file << endl;

        if (op < 0) {
            for (int i = 0; i < max_ch; ++i) {
                int ich = i >= (int) x_data.size () ? x_data.size () - 1 : i;
                int dch = i >= (int) y_data.size () ? y_data.size () - 1 : i;                
                for (int s = 0; s <  out_sz; ++s) {
                    if (s + maxes[i] < (int) x_data[ich].size ()) x_data[ich][s + maxes[i]] -= y_data[dch][s] * mix;
                }
            }

            stringstream res_file;
            res_file << remove_extension ((string) output_file) << ".residual.wav";
            write_wav (res_file.str ().c_str (), x_data, x_header);
            cout << "residual saved to: " << res_file.str () << endl << endl;
        }

        delete [] x_FFT;
        delete [] y_FFT;
        delete [] result_FFT;

        return 0;
    } catch (exception& err) {
        cout << "error: " << err.what () << endl;
    } catch (...) {
        cout << "fatal error; quitting..." << endl;
    }
    return 0;
}

// eof

 