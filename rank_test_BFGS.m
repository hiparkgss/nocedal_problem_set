% let's try 5 by 5 random matrix

dimension = 10;

rng(2017);  % same random number every time
B = rand(dimension); % n by n matrix
B = B' * B;  % making positive definite and symmetric
H = inv(B);

% column vectors
s = rand(dimension, 1);
y = rand(dimension, 1);
rho = 1/(s' * y);

% BFGS with B
a = B * s;
alpha =  - 1/(s' * B * s);
updated_B = alpha * (a * a') + rho * (y * y');
r1 = rank(updated_B);
s1 = sprintf('rank of updated_B is %d', r1);
disp(s1);

% BFGS with H

% inside parenthesis
p1 = eye(dimension) - rho * s * y';
p2 = p1';
updated_H = p1 * H * p2 + rho * (s * s') - H;
r2 = rank(updated_B);
s2 = sprintf('rank of updated_H is %d', r2);
disp(s2);