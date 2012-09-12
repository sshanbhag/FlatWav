% function out = checkfft(S)

NFFT = length(S) / 2;

index1 = 1:NFFT;
index2 = (2*NFFT):-1:(NFFT+1);

out = S(index1) .* S(index2);