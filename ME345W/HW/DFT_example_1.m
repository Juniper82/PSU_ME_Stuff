% Script to create a function containing several frequency components
% and take the DFT and display it.
 close all   %close all open figures
 fs = 1000;              %define sampling frequency in Hz
 N = 1000;               %define number of samples
 t = (0:(1/fs):(N-1)/fs).';   %define time values 
 x1 = 2*sin(2*pi*100*t);  %define tone 1
 x2 = 1*sin(2*pi*250*t);  %define tone 2
 x3 = 0*cos(2*pi*1800*t);  %define tone 3
 x4 = 0*cos(2*pi*2400*t); %define tone 4
 y = 1.5 + x1 + x2 + x3 + x4;    %compound signal (a little DC added as well)
 plot(t,y), grid on, title('Signal Waveform')  %plot the time signal
 Y = fft(y);          %take the discrete Fourier transform of y
 k = (0:length(Y)-1);            %define vector of frequency index values
 figure
 plot(k,abs(Y)), grid on, title('DFT Magnitude')  %plot the magnitudes of the FFT
 xlabel('Frequency Index (k)'),ylabel('|Yk|')
 %plot(k,real(Y)), grid on, title('Real Part of FFT')  %plot the real parts of the FFT
 % figure
 %plot(k,imag(Y)), grid on, title('Imaginary Part of FFT')  %plot the imaginary parts of the FFT
 %  figure
