% let's find the q-1 set of column vectors of s that satisfies the
% condition 

n = 3;
q = 10;
a = rand(n, q-1);

% if the rank of a is 2, then P does not have independent columns.
% a(3, :) = a(1, :) + a(2, :); 

rank(a)
P = [];
for i=1:q-1
    P = [P, h(a(:, i))];
end


det(P)


function output = h(s)
% this function unrolls the upper triangle of rank 1 matrix starting from
% (1,1) >> (1, 2) ... >> (2, 2) >> (2, 3) ...>> (3, 3) ...
% and then add s on top.

% make sure that s is a column vector
s = reshape(s, [], 1);
n = length(s);
r1m = s * s';
container = [];

for i=1:n
    container = [container, r1m(i,i:end)];
end  % end of for

container = reshape(container, [], 1);

output = [s; container];

end


