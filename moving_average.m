function y = moving_average(x, w)
%-----------------------------------------------------------------
% computes sliding moving-average of vector x with window size w
%-----------------------------------------------------------------
k = ones(1, w) ./ w;
% don't want to use normal 0 padding for x vector.  instead, add copies
% of 1st and last values to x vector
% first, force x to row vector
x = x(:)';
% store original length of x
xlen = length(x);
% then, pad vector with w copies of first and last x values
x = [ x(1)*ones(1, w) x x(end)*ones(1, w)];
% convolve
y = conv(x, k, 'same');
% truncate to proper length
y = y(w + (1:xlen));
