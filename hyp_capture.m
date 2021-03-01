function [traj, delta_v] = hyp_capture(goal_id,Planet_mass, Planet_radii, Planet_distances,...
             planet_position,v_arr, pos_arr,arr_time, park_r,v_in,RA,i)
%v_arr è in uscita da sv_from_coe (del pianeta in esame)
%v_in è il vettore velocità che estraiamo da lambert terra2->giove
%pos_arr = orb3 (in questo caso)
%planet position dev'essere all'istante di tempo di entrata della sonda
%origin_coe:
%                h    = angular momentum (km^2/s)
%                e    = eccentricity
%                RA   = right ascension of the ascending
%                       node (rad)
%                incl = inclination of the orbit (rad)
%                w    = argument of perigee (rad)
%                TA   = true anomaly (rad)
%                a    = semimajor axis (km)


sun_masses = 1.989 * 10^30;
G = 6.6742e-20; %[km^3/kg/s^2]

%calcolo la SOI del pianeta 
planet_SOI = (Planet_mass/sun_masses)^(2/5)*Planet_distances;
planet_mu  = G * Planet_mass;
planet_rad = Planet_radii;


rp = planet_rad + park_r;       %perigeo

 

    v_inf =norm(v_arr - v_in);         %v_arr è in uscita da sv_from_coe (del pianeta in esame)
    e = 1 + rp*v_inf^2/planet_mu; % eccentricity
    a = rp/(e-1);                % [km] semi-major axis
    
%raggio target
    delta = rp*sqrt(1+2*planet_mu/(rp*v_inf^2));
%angolo tra desiderato ed effettivo
    half_delta = asin(1/e);     %NOTA: sul libro è 2*asin(1/e)
%velocità iperbole
     v_hyp = sqrt(v_inf^2+2*planet_mu/rp);
% velocity of parking orbit
     vc = sqrt(planet_mu/rp);
%angular momentum
      h = delta*v_inf;
%right ascension of ascending node

    n = sqrt(planet_mu/a^3); % mean motion
    
    %% Trajectory computation
    hyp = zeros(100*24*3600/60,3);
    
    %Angle and direction of the orbit at the entering point
    in_dir = (orbit(1,1:3)-orbit(2,1:3))'; %exit vector: (1,1:3)<-(2,1:3)
    in_angle = deg2rad(atan2d_0_360(in_dir(2),in_dir(1)));

    %Initial point (IL PUNTO DOVE VOGLIO ARRIVARE)
    coe = [h, e, RA, incl, 0, 0];
    [rprova,~] = sv_from_coe(coe, planet_mu); %sv_from_coe mi da posizione e velocità da coe
    alpha = deg2rad(atan2d_0_360(rprova(2),rprova(1)));
    
    %Desired characteristics to align with the interplanetary orbit
    xi_des = in_angle - pi - half_delta;
    alpha_des = pi/2 + xi_des;
    w_des = alpha_des - alpha;
    
        for t=0:60:100*24*3600 
        %Using Kepler's method to compute the point
        M = n*t;
        F = kepler_H(e,M);
        cosf = (cosh(F)-e)/(1-e*cosh(F)); %(Eq. 3.41b) Curtis
        f = acos(cosf);
        

        coe = [h, e, RA, incl, w_des+pi,-f];

        [r,~] = sv_from_coe(coe, planet_mu);
        
        %Sometimes the above algorithm produces NaN elements because of the
        %kepler_H function (it seems not to be able to deal with high
        %numbers)
%         if(any(isnan(r)))
%             if(hyp(2,:) ~= [0 0 0]) %from the third point onward
%                 diff = hyp(counter-1,:)-hyp(counter-2,:);
%                 diff = 60*norm(v_in)*diff/norm(diff);
%                 point = hyp(counter-1,:)' + diff';
%             else %first two points
%                 if (goal_id > 4) %for outer planets/elements
%                     coe = [h, e, RA, incl, w_des+pi,-t/6];
%                 else %for inner planets/elements
%                     coe = [h, e, RA, incl, w_des,t/6];
%                 end
%                 [peri,~] = sv_from_coe(coe, planet_mu);
%                 point = pl_r0' + peri';
%             end
%         else %if kepler_H returned a valid result
%              point = pl_r0' + r';
%         end
%         
    
  
end