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