
r = [earth.radius+200000; 0; 0];      %distanza centro della terra satellite
%r = [1.495978707e11,0,0];
v = [0;sqrt(G*earth.mass/norm(r(1,:),2)) ;0]; %velocità iniziale:    
                                        %ellisse            sqrt(G*M/norm(r(1,:),2))
                                        %iperbole          7.7257e+05    
                                        %parabola        10.9257e+03    
periodo = 2*pi*r(1)/v(2);                                       
timetodeltaV =    5.3100e+03;
                         %7.9650e+03;%periodo*1.5; % [secondi]
                         %6.637e+03 ; %periodo*1.25; % [secondi]
escape_velocity  = sqrt(earth.mu * 2 /r(1));

[orbit1,vel_at_timetodeltav] = timeTrajectory_sec(timetodeltaV, earth.mass, r, v);
plot3(orbit1(1,:)*mtoau+terra(1,1), orbit1(2,:)*mtoau+terra(2,1), orbit1(3,:)*mtoau, 'b' ,'LineWidth', 1);

r = orbit1(:,end);
deltaV = escape_velocity - norm(vel_at_timetodeltav,2);   %[m/s]
v = vel_at_timetodeltav + deltaV*[0;1; 0];

time2 = 500000 % tempo dal istante in cui si effettua il delta v fino all uscita dalla sfera di influenza
[orbit2,vel2] = timeTrajectory_sec(time2, earth.mass, r, v);
plot3(orbit2(1,:)*mtoau+terra(1,1),orbit2(2,:)*mtoau +terra(2,1), orbit2(3,:)*mtoau , 'b' ,'LineWidth', 1);

%                               sfera di influenza della terra
y =linspace(0,6.26,100);
plot3(cos(y)*earth.radius_of_influ + terra(1,1), ...
         sin(y)*earth.radius_of_influ + terra(2,1), ...
         zeros(100), 'r' ,'LineWidth', 1);


r = orbit2(:,end) * mtoau + [terra(1,1); +terra(2,1); 0];
v = vel2* (mtoau/sectodays);
time3 = 40;
[orbit3, finalvel] = longpatch(time3 ,sun.mass, r , v);
plot3(orbit3(1,:), orbit3(2,:), orbit3(3,:) , 'b' ,'LineWidth', 1);

% r = orbit2(:,end) + [terra(1,1)/mtoau; +terra(2,1)/mtoau; 0];
% v = vel2;
% time3 = 4000000;
% [orbit3, finalvel] = timeTrajectory_sec(time3 ,sun.mass, r , v);
% plot3(orbit3(1,:)*mtoau, orbit3(2,:)*mtoau, orbit3(3,:)*mtoau , 'b' ,'LineWidth', 1);
%%               FUNCTIONS

function [distance, finalvel] = longpatch(finalTime,actractorMass,initialPos,initialVel)
    G = 6.67408e-11;            %costante di gravitazione universale
    mtoau = 6.684589e-12;    %fattore di conversione metri unita astronomiche
    sectodays = 0.0000115741;       %fattore di conversione secondi  giorni
    G = (G * mtoau^3) / (sectodays^2);
    r = initialPos;                    %posizione iniziale rispetto all'attrattore
    v = initialVel;
    M = actractorMass;          %massa dell'attrattore
    dt = 1;
    t = 0;
    rplot = zeros(3,finalTime); %prealloching memory for save time
    while t <= finalTime
        rlen = norm(r,2);
        accel = -G*M/rlen^2;
        vnext = v + dt*accel*r/rlen;
        rnext = r + dt*v;
        t = t+dt;
        rplot(:,t) = r;
        v = vnext;
        r = rnext;
    end
    distance = rplot;%*mtoau;
    finalvel = v;%*mtoau;
end

function [distance, finalvel] = timeTrajectory_sec(finalTime,actractorMass,initialPos,initialVel)
    G = 6.67408e-11;            %costante di gravitazione universale
    mtoau = 6.684589e-12;    %fattore di conversione metri unita astronomiche
    r = initialPos;                    %posizione iniziale rispetto all'attrattore
    v = initialVel;
    M = actractorMass;          %massa dell'attrattore
    dt = 1;
    t = 0;
    rplot = zeros(3,finalTime); %prealloching memory for save time
    while t <= finalTime
        rlen = norm(r,2);
        accel = -G*M/rlen^2;
        vnext = v + dt*accel*r/rlen;
        rnext = r + dt*v;
        t = t+dt;
        rplot(:,t) = r;
        v = vnext;
        r = rnext;
    end
    distance = rplot;%*mtoau;
    finalvel = v;%*mtoau;
end

function [distance, finalvel, finaltime] = angleTrajectory(finalAngle,actractorMass,initialPos,initialVel)
    G = 6.67408e-11;            %costante di gravitazione universale
    mtoau = 6.684589e-12;    %fattore di conversione metri unita astronomiche
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
        rplot(:,t) = r;
        v = vnext;
        r = rnext;
        rplot = [rplot, rnext];
        tempDot = dot(initialPos,r);
        costheta = tempDot/ (rlen*rplen);
        theta =real(acos(costheta));
    end
    distance = rplot;%*mtoau;
    finalvel = v;%*mtoau;
    finaltime = t;
end
