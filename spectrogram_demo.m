
fs = 10000;
window = 256;
noverlap = floor(0.95*window)
nfft = 512;


T = 0:(1/fs):2;

X = chirp(T,0,2,fs/2);

[S,F,T,P] = spectrogram(X, window, noverlap, nfft, fs);

surf(T,F,10*log10(P),'edgecolor','none'); axis tight; 

view(0,90);

xlabel('Time (Seconds)'); ylabel('Hz')