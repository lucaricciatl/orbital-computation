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