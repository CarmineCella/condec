import soundfile as sf
import numpy as np
import sys

def next_power_of_2(x):  
    return 1 if x == 0 else 2**(x - 1).bit_length()

def main():
    print("conv.py - convolution based reverberation")

    y, syx = sf.read('pedrotti_dump.wav')
    x, srx = sf.read('sweep.wav')
    # y, syx = sf.read('mix_Bn_Cl_Hn.wav')
    # x, srx = sf.read('ClBb-ord-A3-ff_cut.wav')
    # x, srx = sf.read('Hn-ord-D5-ff.wav')
    # srx = 44100
    # y = np.arange (0, 16)
    # x = np.arange (2, 18)


    scale = 1

    N = next_power_of_2(y.shape[0])   
    ft_y = np.fft.fft (y, N)
    print (ft_y[0:16])
    ft_x = np.fft.fft (x, N)
    print (ft_x[0:16])
    v =  0 # max(abs(ft_x)) * 0.00011
    print ("threshold = ", v)

    print ("deconvolving")
    # resynthesis
    # ir_rebuild = np.fft.ifft(ft_y * np.conj (ft_x) / (v + abs (ft_x) ** 2))
    ft_conv = np.vectorize (complex) (np.zeros ((N,)))
    for i in np.arange (0,N):
        if (abs (ft_x[i]) > v):
            ft_conv[i] = ft_y[i] / ft_x[i]
        else:
            ft_conv[i] = 0
    ir_rebuild = np.fft.ifft(ft_conv) 
    print (np.real (ir_rebuild[0:16]))

    sig_len = y.shape[0] - x.shape[0]
    ir_rebuild = ir_rebuild[1:sig_len] * scale

    print ("saving data")
    sf.write ("ir_rebuild.wav", np.abs (ir_rebuild), srx)

    print ("end")

# main call
main()


