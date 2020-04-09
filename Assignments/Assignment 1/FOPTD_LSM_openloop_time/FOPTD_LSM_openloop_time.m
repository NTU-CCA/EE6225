
%define the original process transfer function
z = [-3.5];
p = [-1,-1,-1,-2,-2,-2,-5];
k = 1;
t_delay = 2.5;
G = zpk(z,p,k,'inputdelay',t_delay);

%define sampleing period which starts from apparent time delay
t_s_start = 4.5;   %sampling start time
t_s_end = 34.5;   %sampling stop time
Ts = 1;      %sampling interval
s_num = (t_s_end-t_s_start)/Ts;   %total sampling number
t_sample = [t_s_start:Ts:t_s_end-Ts];
[y,t] = step(G,t_sample);
 
% calculate psi,gamma for the least squares method
psi = zeros(s_num,3); %initialization
psi(1,1) = Ts/2*y(1);
for i =2:s_num-1
      psi(i,1) = psi(i-1,1) + Ts*y(i);
end
psi(s_num,1) = psi(s_num-1,1) + Ts/2*y(s_num);
psi(:,1) = -psi(:,1);
%trapezoidal integration rule to calculate the intergral of continuous time signal y(¦Ó)
psi(:,2) = -1;%The magnitude of the input step function is 1.
psi(:,3) = 1*t;
gamma = y;

%use LSM to calculate the parameters
theta = ((psi'*psi)^-1)*psi'*gamma;
a1 = theta(1,1);
b1 = theta(3,1);
L  = theta(2,1)/theta(3,1);
 
num = [b1];
den  = [1,a1];
Gn = tf(num,den,'inputdelay',L);

figure('name','FOPTD TIME DOMAIN')

step(G,Gn,35);                        %step test and plot the two responses
figure,nyquist(G,Gn,{0.0001,0.60});   %plot the Nyquist chart in the low frequency

