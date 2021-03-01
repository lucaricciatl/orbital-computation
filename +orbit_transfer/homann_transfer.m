function [distance, finalvel , final_t] = homann_transfer(finalAngle,actractorMass,initialPos,initialVel)
    %%% homan transfer --- final angle [0:pi] %%%

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
    final_t = t 
end



