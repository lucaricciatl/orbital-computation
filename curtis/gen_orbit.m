function [dep_r, dep_v, arr_r, arr_v, flight, orb_oe] = ...
                        gen_orbit(dep_id, arr_id, dep_time, arr_time,flag)
%GENERA L'ORBITA INTERPLANETARIA DELLA NAVICELLA
% GEN_ORBIT(id_pianeta partenza, id_pianeta arrivo, data della partenza, data di arrivo) generates info about
%   the interplanetary orbit of a 1000kg spacecraft departing from 
%   DEP_ID at DEP_TIME and arriving to ARR_ID at ARR_TIME.
%
%   [dep_r,dep_v,arr_r,arr_v,flight,orb_oe] = GEN_ORBIT(...)
%   returns the starting and final position DEP_R, ARR_R
%   of the spacecraft (corresponding to the main bodies position
%   due to the use of the patched conics method), along with its
%   starting and final velocities DEP_V, ARR_V, the time of flight
%   FLIGHT and the orbital elements ORB_OE of the trajectory.
%   
%    dep_id,arr_id - departure/arrival body identifier:
%                1 = Sun
%                2 = Venere
%                3 = Earth
%                4 = Giove
%                5 = Io

%     dep_time     - array specifying time of departure with elements 
%                    (in this order):
%                     year         - range: 1901 - 2099
%                     month        - range: 1 - 12
%                     day          - range: 1 - 31  
%                     hour         - range: 0 - 23  (metteremo 0)
%                     minute       - range: 0 - 60  (metteremo 0)
%                     second       - range: 0 - 60  (metteremo 0)
% 
%     arr_time     - array specifying time of arrival with elements 
%                    (in this order):
%                     year         - range: 1901 - 2099
%                     month        - range: 1 - 12
%                     day          - range: 1 - 31
%                     hour         - range: 0 - 23  (metteremo 0)
%                     minute       - range: 0 - 60  (metteremo 0) 
%                     second       - range: 0 - 60  (metteremo 0)
%
%     orb_oe - orbital elements [h e RA incl w TA a]
%                     where
%                    h    = angular momentum (km^2/s)
%                    e    = eccentricity
%                    RA   = right ascension of the ascending
%                           node (rad)
%                    incl = inclination of the orbit (rad)
%                    w    = argument of perigee (rad)
%                    TA   = true anomaly (rad)
%                    a    = semimajor axis (km)
% 
%     flight       - tempo del volo dal pianeta 1 al 2 (in giorni)
%
%     vinf1, vinf2 - hyperbolic excess velocities at departure
%                    and arrival (km/s)
%
%     dep_r        - position of the main body at departure
%
%     dep_v        - velocity of the spacecraft at departure from
%                    origin SOI
%     arr_r        - position of the main body at arrival
%
%     arr_v        - velocity of the spacecraft at arrival in the 
%                    destination SOI
% 
% User M-functions required: interplanetary, coe_from_sv,
%                             month_planet_names

    %% Argument validation
    validateattributes(dep_time,{'double'},{'size',[1 6]}) % data di partenza
    validateattributes(arr_time,{'double'},{'size',[1 6]}) % data di arrivo

    %% Data
    %Sun mu
    global mu
    global Earth_r0
    global Earth_r1
    global Earth_r2
    global Earth_r3
    global Venere_r0
    global Venere_r1
    global Venere_r2
    global Giove_r2
    global Giove_r3
    deg = pi/180; 

    % Departure pianteta partenza
    planet_id = dep_id;
    year      = dep_time(1);
    month     = dep_time(2);
    day       = dep_time(3);
    hour      = dep_time(4);
    minute    = dep_time(5);
    second    = dep_time(6);
    depart = [planet_id  year  month  day  hour  minute  second];

    % Arrival
    planet_id = arr_id;
    year      = arr_time(1);
    month     = arr_time(2);
    day       = arr_time(3);
    hour      = arr_time(4);
    minute    = arr_time(5);
    second    = arr_time(6);
    arrive = [planet_id  year  month  day  hour  minute  second];

    %% Adding for Venere - Terra 
%--------------------------------------
%                1 = Sun
%                2 = Venere
%                3 = Earth
%                4 = Giove
%                5 = Io

    if (((dep_id == 2 && arr_id == 3)||(dep_id == 3 && arr_id == 4)) && flag ~= 0) %Flag!= 0 fa il flyby
        % Conversions
        au2km = 149597870.700; % [km]
        auday2kms = 149597870.700/86400.0; % [km/s]
         if (flag == 1) %from Venere (o Terra) to added point
            [~, dep_r, dep_pv, ~] = ...
                planet_elements_and_sv(dep_id, ...
                dep_time(1),dep_time(2),dep_time(3),...
                dep_time(4),dep_time(5),dep_time(6));
% modulo della distanza del pianeta venere e velocita' al momento specifico           
            [orb_oe,~,~,~] = planet_elements_and_sv(arr_id, ...
                arr_time(1),arr_time(2),arr_time(3),...
                arr_time(4),arr_time(5),arr_time(6));
% vettore elementi orbitali venere

            % Index calculation based on date
            j0_arr     = J0(arr_time(1),arr_time(2),arr_time(3));
            ut_arr     = (arr_time(4) + arr_time(5)/60 + arr_time(6)/3600)/24;
            jd_arr     = j0_arr + ut_arr;
            j0_dep     = J0(dep_time(1),dep_time(2),dep_time(3));
            ut_dep     = (dep_time(4) + dep_time(5)/60 + dep_time(6)/3600)/24;
            jd_dep     = j0_dep + ut_dep;
            t          = floor((jd_arr - jd_dep +1)/2);%jd - 2454879.5 + 1; % +1 because index starts at 1
% tempo di volo di arrivo all'added point scritto in data giuliana [days]

            if(dep_id == 2 && arr_id == 3)  %venere to earth
                arr_r = Venere_r2;
            elseif (dep_id == 3 && arr_id == 2) % earth to venere
                arr_r = au2km * [-0.5737 0.4473 0.0270];             
            end
            
% down_data restituisce la posizione in funzione del tempo di volo t[days]
            % Setting position, velocity, coe
%             arr_r = au2km * dawn_data(t,1:3);
% arr_r in UA
%           arr_v = auday2kms * dawn_data(t,4:6);
            
            % Time of flight
            flight = t;
            
            % Velocity computation
            [dep_v, arr_v] = lambert(dep_r, arr_r, t*24*3600, 'pro');
            
            pp = 'Dawn';
            
        elseif (flag == 2) %from added point to Earth

             % Data loading
%             dawn_data = dawn_rv();

            % Index calculation based on date
            j0_arr     = J0(arr_time(1),arr_time(2),arr_time(3));
            % Converte la data di partenza in data giuliana (non e' J2000!)
            ut_arr     = (arr_time(4) + arr_time(5)/60 + arr_time(6)/3600)/24;
            % converte ore minuti e secondi in giorni (a noi non interessa)            
            jd_arr     = j0_arr + ut_arr;
            j0_dep     = J0(dep_time(1),dep_time(2),dep_time(3));
            ut_dep     = (dep_time(4) + dep_time(5)/60 + dep_time(6)/3600)/24;
            jd_dep     = j0_dep + ut_dep;
            t          = floor((jd_arr - jd_dep +1)/2); % +1 because index starts at 1
% al tempo t sei nell'added point (pianeta sul quale fai il flyby)

            if(dep_id == 2 && arr_id == 3)  %venere to earth
                dep_r =  Venere_r1;
            elseif (dep_id == 3 && arr_id == 2) % earth to venere
                dep_r = Earth_r3;             
            end
            
%             % Setting position, velocity, coe
%             dep_r = au2km * dawn_data(t,1:3);
% %             dep_v = auday2kms * dawn_data(t,4:6);
            
            % Time of flight
            flight = t;
            
            [orb_oe, arr_r,~,~] = planet_elements_and_sv(arr_id, ...
                arr_time(1),arr_time(2),arr_time(3),...
                arr_time(4),arr_time(5),arr_time(6));
            
            % Velocity computation
            [dep_v, arr_v] = lambert(dep_r, arr_r, t*24*3600, 'pro');
        end
    else
        
        %% Computations
        [planet1, planet2, trajectory] = interplanetary(depart, arrive);
% depart e arrive sono vettori contenenti pianeta di riferimento e data
        %Rp1, Vp1: state vector of planet 1 at departure (km, km/s)
        %R1, V1: heliocentric state vector of the spacecraft 
        %at departure (km, km/s)
        %trajectory vettori velicita' V1 e V2
        R1  = planet1(1,1:3);
        Vp1 = planet1(1,4:6);
        jd1 = planet1(1,7);

        %Rp2, Vp2: state vector of planet 2 at arrival (km, km/s)
        %R2, V2: heliocentric state vector of the spacecraft at
        %arrival (km, km/s)
        R2  = planet2(1,1:3);
        Vp2 = planet2(1,4:6);
        jd2 = planet2(1,7);

        V1  = trajectory(1,1:3);
        V2  = trajectory(1,4:6);

        %time of flight in Julian days
        tof = jd2 - jd1;

        % Use Algorithm 4.2 to find the orbital elements of the
        % spacecraft trajectory based on [Rp1, V1]...
        coe  = coe_from_sv(R1, V1, mu);
        %   ... and [R2, V2]
        coe2 = coe_from_sv(R2, V2, mu);

        % Equations 8.94 and 8.95:
        vinf1 = V1 - Vp1;
        vinf2 = V2 - Vp2;

        
    %% Output info
    % Echo the input data and output the solution to
    % the command window:
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (flag == 1)
        [mm, ~] = month_planet_names(depart(3),depart(1));
        pp = 'Dawn';
    else
        [mm, pp] = month_planet_names(depart(3),depart(1));
    end
        
    fprintf('-----------------------------------------------------')
    fprintf('\n\n Departure:\n');
    fprintf('\n   Planet: %s', pp);%planet_name(depart(1)))
    fprintf('\n   Year  : %g', depart(2))
    fprintf('\n   Month : %s', mm);%month_name(depart(3)))
    fprintf('\n   Day   : %g', depart(4))
    fprintf('\n   Hour  : %g', depart(5))
    fprintf('\n   Minute: %g', depart(6))
    fprintf('\n   Second: %g', depart(7))
    fprintf('\n   Planet position vector (km)    = [%g  %g  %g]', ...
                                                   R1(1),R1(2), R1(3))

    fprintf('\n   Magnitude                      = %g\n', norm(R1))

    fprintf('\n   Planet velocity (km/s)         = [%g  %g  %g]', ...
                                     Vp1(1), Vp1(2), Vp1(3))

    fprintf('\n   Magnitude                      = %g\n', norm(Vp1))

    fprintf('\n   Spacecraft velocity (km/s)     = [%g  %g  %g]', ...
                                                   V1(1), V1(2), V1(3))

    fprintf('\n   Magnitude                      = %g\n', norm(V1))

    fprintf('\n   v-infinity at departure (km/s) = [%g  %g  %g]', ...
                                           vinf1(1), vinf1(2), vinf1(3))

    fprintf('\n   Magnitude                      = %g\n', norm(vinf1))
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fprintf('\n\n  ~~~~~~~~~~~~\n')
    fprintf('\n\n Time of flight = %g days\n', tof)
    fprintf('\n\n  ~~~~~~~~~~~~\n')
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(flag == 2)
        [mm,~] = month_planet_names(arrive(3),arrive(1));
        pp = 'Dawn';
    else
        [mm,pp] = month_planet_names(arrive(3),arrive(1));
    end
    fprintf('\n\n Arrival:\n');
    fprintf('\n   Planet: %s', pp);
    fprintf('\n   Year  : %g', arrive(2))
    fprintf('\n   Month : %s', mm);
    fprintf('\n   Day   : %g', arrive(4))
    fprintf('\n   Hour  : %g', arrive(5))
    fprintf('\n   Minute: %g', arrive(6))
    fprintf('\n   Second: %g', arrive(7))
    fprintf('\n   Planet position vector (km)   = [%g  %g  %g]', ...
                                                  R2(1), R2(2), R2(3))

    fprintf('\n   Magnitude                     = %g\n', norm(R1))

    fprintf('\n   Planet velocity (km/s)        = [%g  %g  %g]', ...
                                      Vp2(1), Vp2(2), Vp2(3))

    fprintf('\n   Magnitude                     = %g\n', norm(Vp2))

    fprintf('\n   Spacecraft Velocity (km/s)    = [%g  %g  %g]', ...
                                                  V2(1), V2(2), V2(3))

    fprintf('\n   Magnitude                     = %g\n', norm(V2))

    fprintf('\n   v-infinity at arrival (km/s)  = [%g  %g  %g]', ...
                                         vinf2(1), vinf2(2), vinf2(3))

    fprintf('\n   Magnitude                     = %g', norm(vinf2))

    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fprintf('\n\n  ~~~~~~~~~~~~\n')
    fprintf('\n\n\n Orbital elements of flight trajectory:\n')

    fprintf('\n  Angular momentum (km^2/s)                   = %g',...
                                                               coe(1))
    fprintf('\n  Eccentricity                                = %g',...
                                                               coe(2))
    fprintf('\n  Right ascension of the ascending node (deg) = %g',...
                                                           coe(3)/deg)
    fprintf('\n  Inclination to the ecliptic (deg)           = %g',...
                                                           coe(4)/deg)
    fprintf('\n  Argument of perihelion (deg)                = %g',...
                                                           coe(5)/deg)
    fprintf('\n  True anomaly at departure (deg)             = %g',...
                                                           coe(6)/deg)
    fprintf('\n  True anomaly at arrival (deg)               = %g\n', ...
                                                          coe2(6)/deg)
    fprintf('\n  Semimajor axis (km)                         = %g',...
                                                               coe(7))
	%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % If the orbit is an ellipse, output the period:
    if coe(2) < 1
        fprintf('\n  Period (days)                               = %g', ...
                                          2*pi/sqrt(mu)*coe(7)^1.5/24/3600)
    end
    fprintf('\n-----------------------------------------------------\n')
    
    %% Output arguments
    dep_r = R1;
    dep_v = V1;
    arr_r = R2;
    arr_v = V2;
    flight = tof;
    %         h     , e     , RA    , incl  , w     , TA    , a
    orb_oe = [coe(1), coe(2), coe(3), coe(4), coe(5), coe(6), coe(7)];
end
end