clear all
close all
clc

G = 6.67408e-11;       %costante di gravitazione universale
M = 5.9722e24;          %massa della terra
m = 1000;                   %massa satellite

r = [6.678e6, 0, 0];      %distanza centro della terra satellite
v = [0,sqrt(G*M/norm(r(1,:),2)),0]; %velocità iniziale
t = 0.0;
dt = 1;
deltaV = 2420;
timetodeltaV = 5431;

rplot = r;

while t <= timetodeltaV
    
    rlen = norm(r,2);
    accel = -G*M/rlen^2;
    vnext = v + dt*accel*r/rlen;
    rnext = r + dt*v;
    t = t + dt;
    
    rplot = [rplot;rnext];
    v = vnext;
    r = rnext;
 
end

hold on
orbit1 = plot(rplot(:,1), rplot(:,2) , 'b' ,'LineWidth', 1);

rplot = r;

theta = 0;
rp = r;
rplen = norm(rp, 2);
vnext = v + deltaV*[0 1 0];
v = vnext;

while theta <= 3.14159
    
    rlen = norm(r,2);
    accel = -G*M/rlen^2;
    
    vnext = v + dt*accel*r/rlen;
    rnext = r + dt*v;
    t = t + dt;
    
    rplot = [rplot; rnext];
    v = vnext;
    r = rnext;
    
    tempDot = dot(rp,r);
    costheta = tempDot/ (rlen*rplen);
    theta = acos(costheta);
    
end

orbit2 = plot(rplot(:,1), rplot(:,2) , 'r' ,'LineWidth', 1);

%%
rplot = r;
tstart = t;

orbitperiod = (2*pi)*sqrt((norm(r,2).^3)/(G*M));
orbitv = sqrt((G*M)/norm(r,2));
deltav = orbitv - norm(v,2);

vnext = v + deltav*[0 -1 0];
v = vnext;

while t< (tstart+orbitperiod)
    
    rlen = norm(r,2);
    accel = -G*M/rlen^2;
    vnext = v + dt*accel*r/rlen;
    rnext = r + dt*v;
    t = t + dt;
    
    rplot = [rplot; rnext];
    v = vnext;
    r = rnext;
    
end


orbit3 = plot(rplot(:,1), rplot(:,2) , 'k' ,'LineWidth', 1);

