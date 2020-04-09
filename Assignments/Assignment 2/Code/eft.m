function[G_hat]=eft(G, Lambda, Gamma )
Am = 3;
G_size = size(G);
for i = 1:G_size(1)
    for j = 1:G_size(2)
    K(i,j) = G(i,j).num{1}(end);
    L(i,j) = G(i,j).InputDelay + G(i,j).ioDelay;
    T(i,j) = G(i,j).den;
    if size(T{i,j}) == [1 2]
        G_hat(i,j)=tf([K(i,j)/Lambda(i,j)],[Gamma(i,j)*T{i,j}(1) 1],'inputdelay',Gamma(i,j)*L(i,j));
    else %size(T{i,j}) == [1 3]
        G_hat(i,j)=tf([K(i,j)/Lambda(i,j)],[T{i,j}(1) T{i,j}(2)*Gamma(i,j) 1],'inputdelay',L(i,j)*Gamma(i,j));  
    end
    end
end
