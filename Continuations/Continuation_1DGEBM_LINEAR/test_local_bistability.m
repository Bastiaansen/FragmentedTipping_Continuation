alpha_1 = 0.7;
alpha_2 = 0.289;
T_1 = 260;
T_2 = 293;
K = 0.1;
A = -339.647;
B = 2.218;

Q0 = 341.3;

%%

Q = @(x) Q0 * (1 - 0.241 * (3*x.^2-1));
alpha = @(T) alpha_1 + (alpha_2-alpha_1) * (1 + tanh(K * (T - (T_1+T_2)/2)))/2;

%%

mu = @(x,T) (A+B*T) - Q(x) .* (1 - alpha(T));

%%

XX = linspace(-1,1,101);
Ts = 200:1:350;

figure()

for i = 1:length(XX)
    plot(mu(XX(i),Ts),Ts)
    hold on
end