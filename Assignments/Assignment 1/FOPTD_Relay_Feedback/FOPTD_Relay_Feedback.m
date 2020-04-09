sim('FOPTD_Relay_Feedback_sim');
plot(simout);

%define the original process transfer function
z = [-3.5];
p = [-1,-1,-1,-2,-2,-2,-5];
k = 1;
t_delay = 2.5;
G = zpk(z,p,k,'inputdelay',t_delay);

%read from the scope output
L = 4.508;
Pu = 13.19;
Kp = 0.0875;
h = 1;
a = 0.07449;

Wu = 2*pi/Pu;
Ku = 4*h/(pi*a);
T = sqrt((Kp*Ku)^2-1)/Wu;

num = [Kp];
den  = [T,1];
Gn = tf(num,den,'inputdelay',L);
figure,step(G,Gn,35);
figure,nyquist(G,Gn,{0.0001,0.56});
