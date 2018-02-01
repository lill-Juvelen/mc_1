%% 2b importance sampling
close all
N=100;
v = linspace(3,25,N);
instgam = gampdf(v,3.1,2.5); % used for plots in document

for  i = 1:length(k)
f(i,:) = wblpdf(v,lambda(i),k(i));
r(i,:) = P(v)'.*(f(i,:)./instgam); % compute the ratio p(v)*(f(v))/g(v)
end

hold on
cmap = colormap(hsv(24));
c = 1;
for  i = 1:length(k)
plot(v,f(i,:), 'LineWidth', 2,'Color', cmap(c, :))
hold on
c = c +3;
end
plot(v,instgam,'Linewidth', 2, 'Color', 'k')
legend('jan & feb', 'march', 'april', 'may', 'june & july & aug', 'sept', 'okt', 'nov & dec', 'g(v)')
title('Probability density functions for each month with instrumental density g (in black)')
xlabel('Windspeed (m/s)')
ylabel('Probability')
%% 
title('The different f_{i}(v) and our chosen instrumental density (*)')
xlabel('Wind speed (m/s)')
legend('jan & feb', 'march', 'april', 'may', 'june & july & aug', 'sept', 'oct', 'nov & dec', 'g(v)')

figure(3)
v = linspace(3,25,N);
cmap = colormap(hsv(24));
c = 1;
for i =1:length(k)
    plot(v, r(i,:), 'Linewidth', 2, 'Color', cmap(c,:))
    hold on
    c = c +3;
end
title('Plot of ratio P(v)f_i(v)/g(v)')
xlabel('Wind speed (m/s)')
legend('feb', 'march', 'april', 'may & july', 'june & aug', 'sept', 'oct', 'nov & dec & jan')

%% Confidence interval

N = 10000;
x = gamrnd(3.1,2.5,N,1);
M = [];
Ilow = [];
Iup = [];
for i = 1:N
    P2 = P(x(1:i)).*wblpdf(x(1:i),lambda(1),k(1))./gampdf(x(1:i),3.1,2.5); %phi*f/g
    m = mean(P2);
    sig = std(P2);
    M = [M;m];
    Ilow = [Ilow;m-lambda(1)*sig/sqrt(i)];
    Iup = [Iup;m+lambda(1)*sig/sqrt(i)];
end

figure;
plot(1:N,M,'r')
hold on
plot(1:N,Ilow,'b')
hold on
plot(1:N,Iup,'b')

%% 2c
N = 1000;
tau = [];
Ilow_AS = [];
Iup_AS = [];

u = rand(1,N/2);
mon = 1; %month

for i = 1:N/2
    V = P(wblinv(u(1:i),lambda(mon),k(mon)));
    Vprim = P(wblinv(1-u(1:i),lambda(mon),k(mon)));
    W = (V+Vprim)./2;
    sig = std(W);
    m = mean(W);
    tau = [tau; m];
    Ilow_AS = [Ilow_AS;m-norminv(0.975)*sig./sqrt(i/2)];
    Iup_AS = [Iup_AS;m+norminv(0.975)*sig./sqrt(i/2)];
end
hold on
plot(1:length(Ilow_AS),Ilow_AS)
plot(1:length(Iup_AS),Iup_AS)
plot(1:length(tau),tau)

%% 2d Calculate probability that the turbine delivers power
%i.e calculate probability that wind is in invterval [3,25]
for i = 1:length(k)
   MonP(i) = integral(@(y)wblpdf(y,lambda(i), k(i)),3,25);
end
bar(1:12,MonP)

%% 2e 
mon = 1;
N = 1000;
rho = 1.225; % air density
d=112; %rotor diameter
Ptot = @(x) rho.*pi./8.*d.^2.*x.^3;

sample = wblrnd(lambda(mon),k(mon),N,1);
rat = [];
Ilow = [];
Iup=[];

for i = 1:N
    ratio = mean(P(x(1:i)))./mean(Ptot(x(1:i)));
    rat = [rat,ratio];
    sig = std(P(x(1:i)))./std(Ptot(x(1:i)));
    Ilow = [Ilow;ratio - lambda(mon)*sig/sqrt(i)];
    Iup = [Iup;ratio + lambda(mon)*sig/sqrt(i)];
end

hold on
plot(1:N,Ilow,'r--');
plot(1:N,Iup,'r--');
plot(1:N,rat,'b');

%% 2f
MonP >= 0.9

mean(MonP)
