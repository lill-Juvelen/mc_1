clear all
close all
% Generates random number from wbl distribution
alpha = 9.2;
beta = 1.7;
rho = 1.225;
d = 90;
sample = wblrnd(alpha,beta,1000,1);

load powercurve_V112
tn = mean(P(sample))

a = 3
b = 25
P_tr = @(x) (P(x) - P(a)) / (P(b) - P(a)) %F(x | x in I)
tn_tr = mean(P(sample))
%P = 1/2 * rho * pi * (d^2)/4 * sample.^3
%tn = mean(P)*0.4 % = 1.8977e+06







figure(1)
subplot(1,2,1)
plot(sample, P(sample), '.')
subplot(1,2,2)
plot(sample, P_tr(sample), '.')


%% jämföra med min instrumental density.. plotta och beräkna ration mellan dem. 
% linspace 0 - 30 och beräkna ration. helst ska insturmental ligga innanför
% x = 3 och x = 25. 
figure(2)
subplot(1,2,1)
hist(sample, 100)
subplot(1,2,2)
isample= wblrnd(14, 2.6, 1000,1)
hist(isample, 100, 'r')

