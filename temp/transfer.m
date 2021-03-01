clear all
close all
clc

%%               MAIN
G = 6.67408e-11;       %costante di gravitazione universale
M = 5.9722e24;          %massa della terra
m = 1000;                   %massa satellite
r = [6.678e6, 0, 0];      %distanza centro della terra satellite
v = [0,sqrt(G*M/norm(r(1,:),2))   ,0]; %velocità iniziale:    
                                        %ellisse            sqrt(G*M/norm(r(1,:),2))
                                        %iperbole          7.7257e+05    
                                        %parabola        10.9257e+03                                  
timetodeltaV = 5431;
%[orbit1,vel_at_timetodeltav] = timeTrajectory(timetodeltaV, M, r, v);

[orbit1,vel_at_timetodeltav] = angleTrajectory(3.14, M, r, v);
hold on
plot(orbit1(:,1), orbit1(:,2) , 'b' ,'LineWidth', 1);


% %  dynamic plot
% for t = 1: timetodeltaV
%     pause(0.0000000000001)
%     plot(orbit1(t,1), orbit1(t,2),'.')
% end

%%               FUNCTIONS

function [distance, finalvel] = timeTrajectory(finalTime,actractorMass,initialPos,initialVel)
    G = 6.67408e-11;            %costante di gravitazione universale
    r = initialPos;                    %posizione iniziale rispetto all'attrattore
    v = initialVel;
    M = actractorMass;          %massa dell'attrattore
    dt = 1;
    t = 0;
    rplot = zeros(finalTime,3); %prealloching memory for save time
    while t <= finalTime
        rlen = norm(r,2);
        accel = -G*M/rlen^2;
        vnext = v + dt*accel*r/rlen;
        rnext = r + dt*v;
        t = t+dt;
        rplot(t,:) = r;
        v = vnext;
        r = rnext;
    end
    distance = rplot;
    finalvel = v;
end

function [distance, finalvel] = angleTrajectory(finalAngle,actractorMass,initialPos,initialVel)
    G = 6.67408e-11;            %costante di gravitazione universale
    r = initialPos;                    %posizione iniziale rispetto all'attrattore
    rplen = norm(initialPos, 2);
    v = initialVel;
    M = actractorMass;          %massa dell'attrattore
    dt = 1;
    t = 0;
    theta = 0;
    
    while theta <= finalAngle
        rlen = norm(r,2);
        accel = -G*M/rlen^2;
        vnext = v + dt*accel*r/rlen;
        rnext = r + dt*v;
        t = t+dt;
        rplot(t,:) = r;
        v = vnext;
        r = rnext;
        rplot = [rplot; rnext];
        tempDot = dot(initialPos,r);
        costheta = tempDot/ (rlen*rplen);
        theta =real(acos(costheta));
    end
    distance = rplot;
    finalvel = v;
end
