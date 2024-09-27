#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <string>
#include <cstdint>

using namespace std;

// Constants
const double PI = 3.14159265358979323846;

// WAV header structure
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

// Function to read WAV file
std::vector<std::vector<double>> read_wav(const char* filename, WAVHeader& header) {
    std::stringstream err;
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        err << "Cannot open WAV file " << filename;
        throw std::runtime_error(err.str());
    }

    // Read the WAV header
    file.read(reinterpret_cast<char*>(&header), sizeof(WAVHeader));

    // Validate WAV format
    if (std::strncmp(header.riff, "RIFF", 4) != 0 ||
        std::strncmp(header.wave, "WAVE", 4) != 0 ||
        std::strncmp(header.fmt, "fmt ", 4) != 0 ||
        std::strncmp(header.data, "data", 4) != 0) {
        err << "Invalid WAV file format.";
        throw std::runtime_error(err.str());
    }

    if (header.audioFormat != 1 && header.audioFormat != 3) { // Only PCM (1) and IEEE float (3) supported
        err << "Unsupported WAV format (not PCM or IEEE float).";
        throw std::runtime_error(err.str());
    }

    // Calculate the number of samples per channel
    int numSamples = header.dataSize / (header.bitsPerSample / 8);
    numSamples /= header.numChannels;
    std::vector<std::vector<double>> channels(header.numChannels, std::vector<double>(numSamples));

    // Read the audio data
    for (int i = 0; i < numSamples; ++i) {
        for (int ch = 0; ch < header.numChannels; ++ch) {
            if (header.bitsPerSample == 16) {
                int16_t sample;
                file.read(reinterpret_cast<char*>(&sample), sizeof(int16_t));
                channels[ch][i] = sample / 32768.0;  // Normalize to [-1.0, 1.0]
            } else if (header.bitsPerSample == 32 && header.audioFormat == 3) {
                float sample;
                file.read(reinterpret_cast<char*>(&sample), sizeof(float));
                channels[ch][i] = sample;  // Already in [-1.0, 1.0]
            } else {
                err << "Unsupported bits per sample: " << header.bitsPerSample;
                throw std::runtime_error(err.str());
            }
        }
    }

    file.close();
    return channels;
}

// Function to write WAV file (not used in feature extraction but included as per requirement)
void write_wav(const char* filename, std::vector<std::vector<double>>& channels, WAVHeader& header) {
    std::stringstream err;
    std::ofstream file(filename, std::ios::binary);

    if (!file) {
        err << "Cannot write WAV file " << filename;
        throw std::runtime_error(err.str());
    }

    // Update header data size
    header.dataSize = 0;
    for (const auto& channel : channels) {
        header.dataSize += channel.size() * (header.bitsPerSample / 8);
    }

    // Write the WAV header
    file.write(reinterpret_cast<char*>(&header), sizeof(WAVHeader));

    int numSamples = channels[0].size();
    for (int i = 0; i < numSamples; ++i) {
        for (size_t ch = 0; ch < channels.size(); ++ch) {
            if (header.bitsPerSample == 16) {
                int16_t intSample = static_cast<int16_t>(std::round(channels[ch][i] * 32767));
                file.write(reinterpret_cast<char*>(&intSample), sizeof(int16_t));
            } else if (header.bitsPerSample == 32 && header.audioFormat == 3) {
                float floatSample = static_cast<float>(channels[ch][i]);
                file.write(reinterpret_cast<char*>(&floatSample), sizeof(float));
            } else {
                err << "Unsupported bits per sample or audio format.";
                throw std::runtime_error(err.str());
            }
        }
    }

    file.close();
}

// Helper function to remove file extension
std::string remove_extension(const std::string& filename) {
    size_t last_dot = filename.find_last_of(".");
    if (last_dot != std::string::npos) {
        return filename.substr(0, last_dot);
    }
    return filename; // No extension found, return the original filename
}

// Function to compute the next power of two
int next_power_of_two(int n) {
    return static_cast<int>(pow(2, ceil(log2(n))));
}

// Recursive FFT function
void fft(std::complex<double>* signal, int N) {
    if (N <= 1) return;

    // Divide
    std::vector<std::complex<double>> even(N / 2);
    std::vector<std::complex<double>> odd(N / 2);
    for (int i = 0; i < N / 2; ++i) {
        even[i] = signal[i * 2];
        odd[i] = signal[i * 2 + 1];
    }

    // Conquer
    fft(even.data(), N / 2);
    fft(odd.data(), N / 2);

    // Combine
    for (int k = 0; k < N / 2; ++k) {
        std::complex<double> t = std::polar(1.0, -2 * PI * k / N) * odd[k];
        signal[k] = even[k] + t;
        signal[k + N / 2] = even[k] - t;
    }
}

// Recursive IFFT function
void ifft(std::complex<double>* signal, int N) {
    // Take conjugate
    for (int i = 0; i < N; ++i) {
        signal[i] = std::conj(signal[i]);
    }

    // Perform FFT
    fft(signal, N);

    // Take conjugate again and divide by N
    for (int i = 0; i < N; ++i) {
        signal[i] = std::conj(signal[i]) / static_cast<double>(N);
    }
}

// Function to generate a Hann window
std::vector<double> hann_window(int size) {
    std::vector<double> window(size);
    for (int i = 0; i < size; ++i) {
        window[i] = 0.5 * (1 - cos(2 * PI * i / (size - 1)));
    }
    return window;
}

// Function to apply a window to the FFT buffer
void apply_window(std::vector<std::complex<double>>& buffer, const std::vector<double>& window) {
    for (size_t i = 0; i < buffer.size(); ++i) {
        buffer[i] *= window[i];
    }
}

// Feature extraction functions

// Compute Energy
double compute_energy(const std::vector<double>& signal) {
    double energy = 0.0;
    for (double sample : signal) {
        energy += sample * sample;
    }
    return energy / signal.size();
}

// Compute Spectral Centroid
double compute_spectral_centroid(const std::vector<double>& spectrum, double sampleRate, int fftSize) {
    double weighted_sum = 0.0, sum = 0.0;
    for (int i = 0; i < fftSize / 2; ++i) {
        weighted_sum += i * spectrum[i];
        sum += spectrum[i];
    }
    return (sum != 0) ? (weighted_sum / sum) * (sampleRate / fftSize) : 0.0;
}

// Compute Spectral Spread
double compute_spectral_spread(const std::vector<double>& spectrum, double centroid, double sampleRate, int fftSize) {
    double variance = 0.0, sum = 0.0;
    for (int i = 0; i < fftSize / 2; ++i) {
        double frequency = i * (sampleRate / fftSize);
        variance += spectrum[i] * pow(frequency - centroid, 2);
        sum += spectrum[i];
    }
    return (sum != 0) ? sqrt(variance / sum) : 0.0;
}

// Compute Spectral Flux
double compute_spectral_flux(const std::vector<double>& spectrum, const std::vector<double>& prev_spectrum) {
    double flux = 0.0;
    for (size_t i = 0; i < spectrum.size(); ++i) {
        flux += pow(spectrum[i] - prev_spectrum[i], 2);
    }
    return sqrt(flux);
}

int main(int argc, char* argv[]) {
    cout << "[features, ver. 0.1]" << endl << endl;
    cout << "audio features extraction" << endl << endl;
    cout << "(c) 2024 Carmine-Emanuele Cella" << endl << endl;

    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <wav_file> <fft_size> <hop_size>\n";
        return 1;
    }

    const char* filename = argv[1];
    int fftSize = std::stoi(argv[2]);
    int hopSize = std::stoi(argv[3]);

    // Validate FFT size (must be power of two)
    if ((fftSize & (fftSize - 1)) != 0) {
        std::cerr << "FFT size must be a power of two.\n";
        return 1;
    }

    try {
        // Read the WAV file
        WAVHeader header;
        std::vector<std::vector<double>> channels = read_wav(filename, header);

        // Mix the channels into a single mono track (mean of channels)
        std::vector<double> mono(channels[0].size(), 0.0);
        for (size_t i = 0; i < mono.size(); ++i) {
            for (int ch = 0; ch < header.numChannels; ++ch) {
                mono[i] += channels[ch][i];
            }
            mono[i] /= header.numChannels;
        }

        // Generate the Hann window
        std::vector<double> window = hann_window(fftSize);

        // Open files for saving features
        std::ofstream energy_file("energy.txt");
        std::ofstream centroid_file("spectral_centroid.txt");
        std::ofstream spread_file("spectral_spread.txt");
        std::ofstream flux_file("spectral_flux.txt");

        if (!energy_file || !centroid_file || !spread_file || !flux_file) {
            std::cerr << "Error opening output files.\n";
            return 1;
        }

        // Initialize previous spectrum for spectral flux calculation
        std::vector<double> prev_spectrum(fftSize / 2, 0.0);

        cout << "proessing..."; cout.flush();
        // Process the audio in frames
        for (size_t start = 0; start + fftSize <= mono.size(); start += hopSize) {
            // Create FFT input buffer with windowing
            std::vector<std::complex<double>> buffer(fftSize);
            for (int i = 0; i < fftSize; ++i) {
                buffer[i] = (start + i < mono.size()) ? mono[start + i] : 0.0;
            }

            // Apply Hann window to the buffer
            apply_window(buffer, window);

            // Perform FFT
            fft(buffer.data(), fftSize);

            // Compute magnitude spectrum
            std::vector<double> magnitude(fftSize / 2, 0.0);
            for (int i = 0; i < fftSize / 2; ++i) {
                magnitude[i] = std::abs(buffer[i]);
            }

            // Compute features
            double energy = compute_energy(magnitude);
            double centroid = compute_spectral_centroid(magnitude, static_cast<double>(header.sampleRate), fftSize);
            double spread = compute_spectral_spread(magnitude, centroid, static_cast<double>(header.sampleRate), fftSize);
            double flux = compute_spectral_flux(magnitude, prev_spectrum);

            // Save features to files
            energy_file << energy << "\n";
            centroid_file << centroid << "\n";
            spread_file << spread << "\n";
            flux_file << flux << "\n";

            // Update previous spectrum
            prev_spectrum = magnitude;
        }

        // Close the output files
        energy_file.close();
        centroid_file.close();
        spread_file.close();
        flux_file.close();

        std::cout << "completed\n";
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
