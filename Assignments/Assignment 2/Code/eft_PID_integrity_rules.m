function[G_hat,PID]=eft_PID_integrity_rules(G, Lambda, Gamma )
Am = 3;
G_size = size(G);
PID = cell(G_size(1),G_size(2));
for i = 1:G_size(1)
    for j = 1:G_size(2)
    K(i,j) = G(i,j).num{1}(end);
    L(i,j) = G(i,j).InputDelay + G(i,j).ioDelay;
    T(i,j) = G(i,j).den;
    if size(T{i,j}) == [1 2]
        if abs(Lambda(i,j)) < 1
            if Gamma(i,j) <= 1 && Gamma(i,j) > 0
                G_hat(i,j)=tf([K(i,j)/Lambda(i,j)],[T{i,j}(1) 1],'inputdelay',L(i,j));
                PID{i,j}={pi*Lambda(i,j)*T{i,j}(1)/2/Am/L(i,j)/K(i,j) pi*Lambda(i,j)/2/Am/L(i,j)/K(i,j) 0};
            else %Gamma(i,j) > 1
                G_hat(i,j)=tf([K(i,j)/Lambda(i,j)],[Gamma(i,j)*T{i,j}(1) 1],'inputdelay',Gamma(i,j)*L(i,j));
                PID{i,j}={pi*Lambda(i,j)*T{i,j}(1)/2/Am/L(i,j)/K(i,j) pi*Lambda(i,j)/2/Am/Gamma(i,j)/L(i,j)/K(i,j) 0};
            end
        else %abs(Lambda(i,j)) >= 1
             if Gamma(i,j) > 1
                G_hat(i,j)=tf([sign(Lambda(i,j))*K(i,j)],[Gamma(i,j)*T{i,j}(1) 1],'inputdelay',Gamma(i,j)*L(i,j));
                PID{i,j}={pi*T{i,j}(1)/2/Am/L(i,j)/K(i,j) pi/2/Am/Gamma(i,j)/L(i,j)/K(i,j) 0};  
             else %Gamma(i,j) <= 1 && Gamma(i,j) > 0
                G_hat(i,j)=tf([sign(Lambda(i,j))*K(i,j)],[T{i,j}(1) 1],'inputdelay',L(i,j));
                PID{i,j}={pi*T{i,j}(1)/2/Am/L(i,j)/K(i,j), pi/2/Am/L(i,j)/K(i,j), 0};
            end
        end
    else %size(T{i,j}) == [1 3]
        if abs(Lambda(i,j)) < 1
            if Gamma(i,j) <= 1
                G_hat(i,j)=tf([K(i,j)/Lambda(i,j)],[T{i,j}(1) T{i,j}(2) 1],'inputdelay',L(i,j));
                PID{i,j}={pi*Lambda(i,j)*T{i,j}(2)/2/Am/L(i,j)/K(i,j) pi*Lambda(i,j)/2/Am/L(i,j)/K(i,j) pi*Lambda(i,j)*T{i,j}(1)/2/Am/L(i,j)/K(i,j)};
            else % Gamma(i,j) > 1
                G_hat(i,j)=tf([K(i,j)/Lambda(i,j)],[T{i,j}(1) T{i,j}(2)/Gamma(i,j) 1],'inputdelay',L(i,j)*Gamma(i,j));
                PID{i,j}={pi*Lambda(i,j)*T{i,j}(2)/2/Am/Gamma(i,j)^2/L(i,j)/K(i,j) pi*Lambda(i,j)/2/Am/Gamma(i,j)/L(i,j)/K(i,j) pi*Lambda(i,j)*T{i,j}(1)/2/Am/Gamma(i,j)/L(i,j)/K(i,j)};
            end
        else %abs(Lambda(i,j)) >= 1
             if Gamma(i,j) > 1
                G_hat(i,j)=tf([sign(Lambda(i,j))*K(i,j)],[T{i,j}(1) T{i,j}(2)/Gamma(i,j) 1],'inputdelay',Gamma(i,j)*L(i,j));
                PID{i,j}={sign(Lambda(i,j))*pi*T{i,j}(2)/2/Am/Gamma(i,j)^2/L(i,j)/K(i,j) sign(Lambda(i,j))*pi/2/Am/Gamma(i,j)/L(i,j)/K(i,j) sign(Lambda(i,j))*pi*T{i,j}(1)/2/Am/Gamma(i,j)/L(i,j)/K(i,j)};
             else % Gamma(i,j) <= 1 
                G_hat(i,j)=tf([sign(Lambda(i,j))*K(i,j)],[T{i,j}(1) T{i,j}(2) 1],'inputdelay',L(i,j));
                PID{i,j}={pi*sign(Lambda(i,j))*T{i,j}(2)/2/Am/L(i,j)/K(i,j) pi*sign(Lambda(i,j))/2/Am/L(i,j)/K(i,j) pi*sign(Lambda(i,j))*T{i,j}(1)/2/Am/L(i,j)/K(i,j)};
            end
        end
    end
    end
end
