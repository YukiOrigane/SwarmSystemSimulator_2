
L = [1 -1 0 0 0 0 0; -1 2 -1 0 0 0 0; 0 -1 2 -1 0 0 0; 0 0 -1 2 -1 0 0; 0 0 0 -1 2 -1 0; 0 0 0 0 -1 2 -1; 0 0 0 0 0 -1 1];
[P,D] = eig(L)

Dinv = pinv(D)

L_hat = P*Dinv*P.'

input = [1;0;0;0;0;0;0]

answer = L_hat*input

plot(answer,'*')

f = fit((1:7).', answer, 'poly2')

hold on
plot(f)