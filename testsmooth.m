

a = 2 * sin(2 * pi * 2 * [0:0.025:1]);
a = a + 2 * (rand(size(a)) -1);

SmoothVal1 = 3;
SmoothVal2 = 5;

figure(2)
subplot(211)
a_sm = conv(a, ones(1, SmoothVal1) ./ SmoothVal1, 'same');
plot([a' a_sm'], '.-')

subplot(212)
b_sm = sgolayfilt(a, SmoothVal1, SmoothVal2);
plot([a' b_sm'], '.-')

