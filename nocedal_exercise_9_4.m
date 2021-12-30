% let's suppose the line is y=2x+1
x1 = (1:6)';
x2 = 2 * x1 + 1;
x0 = ones(6, 1);
x1_sq = x1 .* x1;
x2_sq = x2 .* x2;
x1x2 = x1 .* x2;

P = [x0, x1, x2, x1_sq, x2_sq, x1x2];
rref(P)


% part b
theta = pi/4 * (1:6)';
y1 = cos(theta);
y2 = sin(theta);
y0 = ones(6, 1);
y1_sq = y1 .* y1;
y2_sq = y2 .* y2;
y1y2 = y1 .* y2;

Q = [y0, y1, y2, y1_sq, y2_sq, y1y2];
rref(Q)
