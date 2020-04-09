function[PID]=PID_controller(G)
Am = 3;
G_size = size(G);
PID = cell(G_size(1),G_size(2));
for i = 1:G_size(1)
    K(i,i) = G(i,i).num{1}(end);
    L(i,i) = G(i,i).InputDelay + G(i,i).ioDelay;
    T(i,i) = G(i,i).den;
    if size(T{i,i}) == [1 2]
        PID{i,i}={pi*T{i,i}(1)/2/Am/L(i,i)/K(i,i), pi/2/Am/L(i,i)/K(i,i), 0};
    else %size(T{i,i}) == [1 3]
        PID{i,i}={pi*T{i,i}(2)/2/Am/L(i,i)/K(i,i) pi/2/Am/L(i,i)/K(i,i) pi*T{i,i}(1)/2/Am/L(i,i)/K(i,i)};
    end
    end
end
