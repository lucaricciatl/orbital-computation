function [a, e, i, omega, w, theta] = pv_to_op(R,V,mu)
    %% Calculate all six orbital elements fromega the initial position and velocity
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
    
    clear r v vr H h i k N n E e omega w theta; clc;
    
    %% Set up the initial conditions
    if nargin == 2
        mu = 132712440018;
    end
    
    r = norm(R);
    v = norm(V);
    vr = dot(R,V) / r; % If vr > 0 object is moving away fromega perigee
    
    %% Calculate the specific angular momegaentum
    
    H = cross(R,V)  %momento angolare dell'orbita di trasferimento che si 
                    %trova sull'asse z del sistema solidale all'orbita di trasferimento
    
    h = norm(H)       
    
    %% Calculate node line vector
    
    K = [0, 0, 1]  %versore lungo z
    J = [0 1 0] % versore lungo y
    I = cross(K,H)        %
    N = cross(K,H) / norm(cross(K,H))   %sarebbe n cappuccio
    
    %% Calculate the eccentricity


    E = ((v^2)  - (mu / r) ) %energy
    e = sqrt(1 + ((2 * (E) * h^2) / mu^2)) ;

    %E = (1 / mu)*((((v)^2) - (mu / r)* R) - vr * V)  %energia orbitale specifica
    %e = norm(E)  %
    
    %% Calculate the inclination
    
    i = acos(dot(K,H)/h) ;
    
    
    %% Calculate the right ascension of the ascending node
    k_ = [0, 0, 1];
    N_ = cross(k_,H);
    n_ = norm(N_);
    if N_(2) >= 0
        omega = acos(N(1) / n_) ;
    else
        omega = 2*pi - acos(N_(1) / n_);
    end
    
    
%     if dot(N,J) >= 0
%         omega = acos(dot(I,N)) ;
%     else
%         omega = 2*pi - acos(dot(I,N)) ;
%     end
    
    
    
    %% Calculate the argument of perigee(w)
    k = [0, 0, 1];
    E = (1 / mu)*((v^2 - (mu / r))*R - vr*V);
    e = norm(E);
    N = cross(k,H);
    n = norm(N);
    if E(3) >= 0
        w = acos(dot(I*e)/(e));
    else
        w = 2*pi - acos(dot(N,E)/(n*e)) ;
    end
    
    %% eccentricity vector
    omega 
    i
    w
    Tt = T_matrix.Tzxz(omega,i,w);
    disp(Tt(1:3,1:3));
    x_local = [1 ; 0 ;0];
    disp(Tt(1:3,1:3))
    Tinv = (Tt^-1)
    T = Tinv(1:3,1:3)
    e_vect = e * T*x_local
    
    %% Calculate the true anomegaaly
    
    if vr >= 0
        theta = acos(dot(E/e,R/r)) ;
    else
        theta = 2*pi - acos(dot(E/e,R/r)) ;
    end
    
    E = ((v^2)/2  - (132712440018 / r) )
    a = 132712440018/(2*E); %[km]
end