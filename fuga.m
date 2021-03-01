function [orbit, deltav , days_to_exit, r3, v3]   =   fuga(massa, raggio_influenza, raggio_pianeta, altitudine, orbita_attrattore)
4%
%INPUT
% orbita ->  matrice 3 x sec_to_exit contenente lo storico delle posizioni in terna attrattore
%deltav -> deltav minimo necessario per uscire dalla sfera di influenza
%days_to_exit -> tempo necessario per uscire dalla sfera di influenza
%r3 -> ultima posizione della sonda in terna sole
%v3 -> vettore velocità della sonda alla posizione r3 considerata come
%         somma della velocita relativa al pianeta più la velocità di trascinamento
%         con esso
%OUTPUT
%massa -> massa dell'attrattore
%raggio_influenza, 
%raggio_pianeta, 
%altitudine, 
%orbita_attrattore -> matrice 3 x numeroGiorni contenente lo storico delle posizioni del
%                           pianeta in terna sole.
%                           le manovre di fuga durano solitamente qualche giorno

G = 6.67408e-11;                      %costante di gravitazione universale  
mtoau = 6.684589e-12;              %fattore di conversione metri unita astronomiche
sectodays = 0.0000115741;       %fattore di conversione secondi  giorni
sun_mass = 1.98855e30;
mu = massa*G;


% orbita geocentrica della sonda
r1 = [raggio_pianeta+altitudine; 0; 0];                                 %distanza dal centro della terra al satellite
v1 = [0; sqrt(G*massa/norm(r1,2)) ;0];                  %velocità iniziale
periodo = round(2*pi*norm(r1,2)/norm(v1,2));               %periodo di rivoluzione nell'orbita di parcheggio
time1 = round(periodo*1.5);
[orbit1,v11] = timeTrajectory_sec(time1,massa, r1, v1);

% iperbole di fuga
r2 = orbit1(:,end);
escape_velocity  = sqrt(mu * 2 /norm(r1,2));      %modulo velocità di fuga
v2 = escape_velocity*v11/(norm(v11,2));
deltav =norm(v2,2)-norm(v11);
[orbit2, v21 , time2] = exit_from_soi(raggio_influenza, massa,r2,v2);

% calcolo velocità e posizione della sonda rispetto al sole nell'istante in
% cui essa esce dalla sfera di influenza della terra
days_to_exit = round((time1 + time2 ) /(60*60*24));
earth_at_exit = orbita_attrattore(:,days_to_exit)/mtoau;
r_hat = earth_at_exit/norm(earth_at_exit,2);    %versore posizione terra
alpha = atan2(r_hat(2),r_hat(1))+pi/2;
j = [cos(alpha);sin(alpha);0];
v_terra = sqrt(G*sun_mass/norm(earth_at_exit,2))*j ;   %vettore velocità della terra

r3 = (orbit2(:,end) + earth_at_exit) * mtoau;
v3 =( v21 + v_terra ) * mtoau/sectodays;           %velocità relativa alla terra + velocità di trascinamento con la terra

orbit = [orbit1,orbit2];

figure 
hold on
   axis tight         %
     set(gca,'Color','k')
    set(gcf,'Color','black');
    set(gca,'xcolor','w')
    set(gca,'ycolor','w')
    set(gca,'zcolor','w')
    grid on
    daspect([1 1 1])        %
    axis([-10*10^4 10*10^4 -10*10^4 10*10^4 -10*10^4 10*10^4])
% sfera di influenza della terra
%     y =linspace(0,6.26,100);
%     plot3(cos(y)*raggio_influenza/1000, sin(y)*raggio_influenza/1000,zeros(100)/1000, 'r' ,'LineWidth', 1);
    
% orbita geocentrica della sonda

        %% Figure plot
        obj_id = 2;
        body_pos = [0,0,0];
        body_sphere(obj_id,body_pos);
plot3(orbit1(1,:)/1000, orbit1(2,:)/1000,orbit1(3,:)/1000, 'w' ,'LineWidth', 1);
% iperbole di fuga
plot3(orbit2(1,:)/1000, orbit2(2,:)/1000,orbit2(3,:)/1000, 'w' ,'LineWidth', 1);
hold off

return

%%               FUNCTIONS
function [distance, finalvel , finaltime] = exit_from_soi(actractorROI, actractorMass,initialPos,initialVel)
    time = 0;
    pos = initialPos;
    vel = initialVel;
    while norm(actractorROI,2)>norm(pos(:,end),2)
        time = time+1;
        [pos,vel] = timeTrajectory_onlyLastPos(3600, actractorMass, pos, vel);
    end
    finaltime = time*3600;
    [totalpos,lastvel] = timeTrajectory_sec(finaltime, actractorMass, initialPos, initialVel);
    distance = totalpos;
    finalvel = lastvel;
end

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

function [distance, finalvel] = timeTrajectory_onlyLastPos(finalTime,actractorMass,initialPos,initialVel)
    G = 6.67408e-11;            %costante di gravitazione universale
    r = initialPos;                    %posizione iniziale rispetto all'attrattore
    v = initialVel;
    M = actractorMass;          %massa dell'attrattore
    dt = 1;
    t = 0;
    while t <= finalTime
        rlen = norm(r,2);
        accel = -G*M/rlen^2;
        vnext = v + dt*accel*r/rlen;
        rnext = r + dt*v;
        t = t+dt;
        v = vnext;
        r = rnext;
    end
    distance = r;%*mtoau;
    finalvel = v;%*mtoau;
end

function [distance, finalvel] = timeTrajectory_sec(finalTime,actractorMass,initialPos,initialVel)
    G = 6.67408e-11;            %costante di gravitazione universale
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