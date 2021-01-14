%define the original process transfer function
z = [-3.5];
p = [-1,-1,-1,-2,-2,-2,-5];
k = 1;
t_delay = 2.5;
G = zpk(z,p,k,'inputdelay',t_delay);

%define sampleing period which starts from apparent time delay
t_s_start = 0;   %sampling start time
t_s_end = 20;   %sampling stop time
Ts = 0.5;      %sampling interval
s_num = (t_s_end-t_s_start)/Ts;   %total sampling number
t_sample = [t_s_start:Ts:t_s_end-Ts];
[y,t] = step(G,t_sample);

%trapezoidal integration rule
syms w;

% trapz_sin_part = Ts/2*(y(1)-y(s_num))*sin(w*t(1));
% for n =2:s_num-2
%     trapz_sin_part = trapz_sin_part + Ts*(y(n)-y(s_num))*sin(w*t(n));
% end
% trapz_sin_part = trapz_sin_part + Ts/2*(y(s_num-1)-y(s_num))*sin(w*t(s_num-1));
% 
% trapz_cos_part = Ts/2*(y(1)-y(s_num))*cos(w*t(1));
% for n =2:s_num-2
%     trapz_cos_part = trapz_cos_part + Ts*(y(n)-y(s_num))*cos(w*t(n));
% end
% trapz_cos_part = trapz_cos_part + Ts/2*(y(s_num-1)-y(s_num))*cos(w*t(s_num-1));
% 
%  g(w) = (y(s_num)+w*trapz_sin_part)+1i*w*trapz_cos_part;
 

g(w) = (y(s_num)+w*trapz(t,(step(G,t)-y(s_num)).*sin(w*t)))+1i*w*trapz(t,(step(G,t)-y(s_num)).*cos(w*t));

%initialization
%number of the frequency response to be identified 
%in the frequency range (0, ¦Ø c ) 
M = 10;  
fai = zeros(1,M);
w = zeros(1,M);
psi = zeros(M,2);
gamma = zeros(M,1);

W(1) = 0.0; W(2) = 0.001;
fai(1) = 0.0; fai(2) = angle(g(W(2)));

%recursive solution
%calculate the argument and counterpart frequency of the sampling
for n=3:M
    W(n) = W(n-1)-((n-1)*pi/(M-1)+fai(n-1))*(W(n-1)-W(n-2))/(fai(n-1)-fai(n-2));
    fai(n) = phase(g(W(n)));
end

%calculate parameters a1 and b1 by LSM method.
for n=1:M
    psi(n,1) = -((real(g(W(n))))^2+(imag(g(W(n))))^2);
    psi(:,2) = 1;
    gamma(n) = (W(n)^2)*((real(g(W(n))))^2+(imag(g(W(n))))^2);
end
theta = ((psi'*psi)^-1)*psi'*gamma;
a1 = sqrt(theta(1));
b1 = sqrt(theta(2));

%calculate parameters L by LSM method.
psi2 = W';
gamma2 = zeros(M,1);
for n=1:M
    gamma2(n) = -fai(n)-atan(W(n)/a1);
end
L = ((psi2'*psi2)^-1)*psi2'*gamma2;

num = [b1];
den  = [1,a1];
Gn = tf(num,den,'inputdelay',L);
step(G,Gn,35)                    %do the step test and plot
figure, nyquist(G,Gn,{0.0001,0.56});        %plot the Nyquist chart and compare
