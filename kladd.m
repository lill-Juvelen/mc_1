%% 2a)  Create two different confidence intervals for each months
clear all; close all; load powercurve_V112
lambda = [9.7, 9.2, 8, 7.8, 8.1, 9.1, 9.9, 10.6]
k = [2, 2, 1.9, 1.9, 1.9, 2, 1.9, 2]
N = 10000;
kvantil = norminv(1 - 0.05/2)

rng(36)
% Creates confidence interval for produced power for each month
% using crude Monte Carlo method
for i = 1:length(k)
    sample = wblrnd(lambda(i),k(i),N,1); % the random samples
    tn(i) = mean(P(sample)) % expected value estimated from our random samples
    s(i) = std(P(sample)) % standard deviation estimated from our random samples
    % confidence intervalls for all months
    CI_1(i, :) = [tn(i) - s(i)/sqrt(N) * kvantil, tn(i) + s(i)/sqrt(N) *kvantil] 
end

% Creates confidence interval for produced value P(v) using samples generated by
% f(v) conditioned on the interval I = (3,25)

% generate random numbers from a conditioned, truncated distributions
% and create a CI for each month

for i = 1:length(k)
    A(i) = wblcdf(3,lambda(i), k(i)); % F(a)
    B(i) = wblcdf(25,lambda(i), k(i)); % F(b)
    Fhat_inv = @(u) wblinv(u*(B(i) - A(i)) + A(i), lambda(i), k(i));
    f_sample = inverseMethod(Fhat_inv, N); % f_sample is from conditioned 
    % distribution. Now we have x from the correct interval but need to
    % compensate the probability in order to actually consider the correct
    % distribution of wind speeds. 
    tn_tr(i) = mean(P(f_sample))*(B(i)-A(i)); % new mean
    s_tr(i) = std(P(f_sample)); % new standard deviation
    CI_tr(i,:) = [tn_tr(i) - kvantil*s_tr(i)/sqrt(N) , tn_tr(i) + kvantil* s_tr(i)/sqrt(N)]; %
end

%% Comparing methods plot
% comparing the two different methods vs numbers of samples for june
% (mon = 5)
mon = 5;
Ilow = []
Iup = []
Ilow_prim = []
Iup_prim = []
tic;
for i = 200:N
    tau = mean(P(sample(1:N)));
    tau_prim = mean(P(f_sample(1:N)))*(B(mon) - A(mon));
    s = std(P(sample(1:N)));
    s_prim = std(P(f_sample(1:N)));
    Ilow(i) = tau - kvantil*s/sqrt(i);
    Ilow_prim(i) = tau_prim - kvantil*s_prim/sqrt(i);
    Iup(i) = tau + kvantil*s/sqrt(i);
    Iup_prim(i)  = tau_prim + kvantil*s_prim/sqrt(i);
    t(i) = tau;
    tprim(i) = tau_prim;
end
toc;

figure(1)
hold on
plot(1:N,Ilow,'r--');
hold on
plot(1:N,Iup,'r--');
hold on
plot(1:N,t,'b');
hold on
plot(1:N,Ilow_prim,'g--');
hold on
plot(1:N,Iup_prim,'g--');
hold on
plot(1:N, tprim, 'k');
hold on


%% histogram plots of f(v) and f(v) conditioned on I = (3, 25). 
figure(5) 
subplot(1,2,1)
hist(sample, 100)
title('Histogram samples from f(v), december')
subplot(1,2,2)
hist(f_sample, 100)
title('Histogram samples from f(v) conditioned on the interval I = (3,25) for december')