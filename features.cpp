#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <complex>
#include <stdexcept>
#include <string>

#include "dsp.h"
#include "utils.h"

using namespace std;

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
