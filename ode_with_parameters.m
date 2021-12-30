x0 = [1; 0];
tspan = [0, 6];
options = odeset('stat','on');
parameters.freq = 20;


[t, x] = ode15s(@(t, x) myode(t, x, parameters), tspan, x0, options);

plot(t, x);

function dxdt = myode(t, x, p)
freq = p.freq;
dxdt = [-freq * x(2); freq * x(1)];

end


