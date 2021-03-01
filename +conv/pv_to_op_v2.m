function [a, e, i, omega, w, theta] = pv_to_op(R,V,mu)
    % Calculate all six orbital elements fromega the initial position and velocity
    %
    %
    % Jeremy Penn
    % 15 October 2017
    %
    % Revision: 15/10/17
    %
    %           20/10/17 - Changed notation to conform with "Orbital
    %           Mechanics for Engineering Students" by Howard D. Curtis.
    %           
    % function orbitElements(R,V,mu)
    %
    % Purpose:  This function calculates the classic orbital elements.
    % 
    % Inputs:   o R  - A 1x3 vector describing the initial position of the
    %                  satellite.[km]
    %           o V  - A 1x3 vector describing the initial velocity of the
    %                  satellite.[km/s]
    %           o mu - Standard gravitationl parameter of the central body
    %                  [OPTIONAL]. Defaults to Earth (398600 [km^3/s^2])
    %
    % Output:   o h     - Specific angular momentum
    %           o e     - eccentricity
    %           o i     - orbital inclination
    %           o omega - right ascension of the ascending node
    %           o w     - argument of perigee
    %           o theta - true anomaly
    %
    
    % Set up the initial conditions
    
    r = norm(R,2);
    v = norm(V,2);
    vr = dot(R,V) / r; % If vr > 0 object is moving away fromega perigee
    
    % Calculate the specific angular momegaentum
    
    H_hat = cross(R,V)  %momento angolare dell'orbita di trasferimento che si 
                            %trova sull'asse z del sistema solidale all'orbita di trasferimento
    h = norm(H,2)       
    
    % Calculate node line vector
    
    K = [0, 0, 1];  %versore lungo z
    J = [0, 1, 0]; % versore lungo y
    I = cross(K,H_hat);        %
    N_hat  = cross(K,H_hat) / norm(cross(K,H_hat));   
    
    % Calculate the eccentricity
    E = ((norm(v,2)^2)/2  - (mu / norm(r,2)) ); %energy
    e = sqrt(1 + ((2 * (E) * h^2) / mu^2));
    
    % Calculate the inclination
    i = acos(dot(K,H_hat)/h) * (180/pi)
    
    % Calculate the right ascension of the ascending node
    if dot(N_hat,J) >= 0
        omega = acos(dot(I,N_hat)) * (180/pi);
    else
        omega = 360 - acos(dot(I,N_hat)) * (180/pi);
    end
    
    %  Calculate the eccentricity vector
    y_local = Rz(omega*(pi/180)+pi)*J;
    e_vect = e*cross(H_hat, y_local)/h
    
    % Calculate the argument of perigee
    if cross(E,K) >= 0
        w = acos(cross(N_hat, e_vect)/e) * (180/pi);
    else
        w = 360 - acos(cross(N_hat, e_vect)/e) * (180/pi);
    end
    
    % Calculate the true anomegaaly
    
    if vr >= 0
        theta = acos(cross(e_vect,R)/e*r) * (180/pi);
    else
        theta = 360 - acos(cross(e_vect,R)/e*r) * (180/pi);
    end
    a = -mu/(2*E); %[km]
end

