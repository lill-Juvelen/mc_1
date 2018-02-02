%% 2c Antithetic sampling
% provar för mars och den blir lite för högt väntevärde?
N = 22000;
tau = [];
Ilow_AS = [];
Iup_AS = [];

u = rand(1,N/2);
mon = 2; %month

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
figure(1)
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
