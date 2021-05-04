# Dynamic-Mechanical-Analysis

Python code to analyse data generated from Optical Interoferometry-Based Micropipette Aspiration. This method enables accurate identification of mechanical properties of materials. The data is generated from a DMA-based Labview model and is saved as TDMS file. In this DMA, pressure is applied in a sum of multiple frequencies.

In order to filter out the single-frequency signal, the maxima of the FFT (def_FFT) are used to find the proper frequencies and a lock-in amplifier principle to find the corresponding amplitudes and phases (def_lock_in).

Multiple_DMA is the executor of the analysis. The output represents the input signal compared to the reference wave based on the found frequencies, amplitudes and phases to see if these variables match with the original signal. 

Eventually the frequency-dependent viscoelastic response of the material will be determined (not yet implemented in the code).



