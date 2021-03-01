function [orbit, total_time, last_vel, total_deltaV]= homan(massa, raggio_pianeta, altitudine)

G = 6.67408e-11;                      %costante di gravitazione universale  
mtoau = 6.684589e-12;                 %fattore di conversione metri unita astronomiche
sectodays = 0.0000115741;             %fattore di conversione secondi  giorni

r1 = [raggio_pianeta+altitudine; 0; 0];      %distanza centro della terra satellite
v1 = [0; sqrt(G*massa/norm(r1,2)) ;0];                  %velocità iniziale
periodo = round(2*pi*norm(r1,2)/norm(v1,2));               %periodo di rivoluzione nell'orbita di parcheggio
% orbita piccola 
[orbit1,v11] = timeTrajectory_sec(periodo, massa, r1, v1);

% trasferimento nella seconda orbita per mezzo del primo deltav
r2 = orbit1(:,end);
deltav1 = 831;
v2 = v11 -deltav1*[0 ;1; 0];
[orbit2,v21, time] = angleTrajectory(3.14, massa, r2, v2);

r3 = orbit2(:,end);
orbitperiod = round((2*pi)*sqrt((norm(r3,2).^3)/(G*massa)));  %9.0370e+04  ->  9e+04
orbitv = sqrt((G*massa)/norm(r3,2));
deltav2 = orbitv - norm(v21,2);
v3 = v21 + deltav2*[0; -1; 0];
[orbit3,v31] = timeTrajectory_sec( orbitperiod, massa, r3, v3);

total_time = (orbitperiod + periodo + time)*sectodays;
orbit_m = [orbit1, orbit2, orbit3];
orbit = orbit_m/1000;

total_deltaV = deltav1 + deltav2;
last_vel = v31;
%-0.2627
% plot3(orbit1(1,:), orbit1(2,:), orbit1(3,:), 'b' ,'LineWidth', 1);
% plot3(orbit2(1,:), orbit2(2,:), orbit2(3,:) , 'b' ,'LineWidth', 1);
% plot3(orbit3(1,:), orbit3(2,:), orbit3(3,:), 'b' ,'LineWidth', 1);
plot3(orbit(1,:), orbit(2,:), orbit(3,:), 'b' ,'LineWidth', 1);
hold off

return

%%               FUNCTIONS


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

end