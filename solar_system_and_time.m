%% traiettorie nel tempo

param
days = 951;                 % durata della simuazione in giorni
time_at_j2000 = 2451545;
days_in_a_century = 36525;

% preallocation for save time
terra = zeros(3,days+1);
venere = zeros(3,days+1);

giorni = 10842; % giorni trascorsi dall'1/1/2000  al 7/9/2029(data di partenza)
initial_time = time_at_j2000 + giorni;
final_time =  initial_time + days;

i=0;
for Teph = initial_time : final_time 
    i = i+1;
    T = (Teph - time_at_j2000) /  days_in_a_century; %  [century]
    
    terra(:,i) = advance(earth, T);
    venere(:,i) = advance(venus,T);
end

[orbit1, total_time, last_vel, total_deltaV]= homan(earth.mass, earth.radius, 200000);
[orbit, deltav , days_to_exit, r3, v3]   =   fuga(earth.mass, earth.radius_of_influ, earth.radius, 200000, terra);
             

% figure
% title('planet and satellite orbit')
% hold on
% plot3(0,0,0,'-*')
%     plot3(terra(1,1), terra(2,1),terra(3,1),'-o')
%     plot3(venere(1,373), venere(2,373),venere(3,373),'-s')
%     plot3(terra(1,end), terra(2,end),terra(3,end),'-*')
%     plot3(terra(1,:), terra(2,:),terra(3,:), 'k' ,'LineWidth', 1)
%     plot3(venere(1,:), venere(2,:),venere(3,:), 'k' ,'LineWidth', 1)
%     plot3(orbit3(1,:), orbit3(2,:),zeros(1,301),'b','LineWidth', 1)
% hold off

% figure 
% hold on
% plot3(-1.5,-1.5,0,'.')
% plot3(1.5,1.5,0,'.')
% for t = 1: days+1
%     pause(0.0001)
%     plot3(terra(1,t), terra(2,t),terra(3,t), '.')
%     plot3(venere(1,t), venere(2,t),venere(3,t),'.')
% end
% hold off


%%                      FUNCTIONS

function pos = advance(planet, T)

    om_bar = planet.om_bar ;%+ planet.om_bar_dot*T;
    om_big = planet.om_big ;%+ planet.om_big_dot*T;
    L = planet.L + planet.L_dot*T;
    e = planet.e ;%+ planet.e_dot*T;
    I =  planet.I ;%+ planet.I_dot*T;
    a = planet.a ;%+ planet.a_dot*T ;
    
    om = om_bar - om_big;               %periheion
    M = L - om_bar;                         %mean anomaly
    
    E =  eccentric_anomaly(M, e);
    
    position_first = [ a*(cos(E)-e);
                             a*sqrt(1-e^2)*sin(E);
                             0 ];
                                    
    Rzxz = Rz( -om_big )*Rx(-I)*Rz(-om);
    
    pos =   Rzxz *position_first; %position_first;
end

function E = eccentric_anomaly(M, e)

    E = M + e*sin(M);
    delta_E = 1;
    tol = 1e-06;
    
    while(abs(delta_E)>tol)
        delta_M = M-(E-e*sin(E));
        delta_E = delta_M/(1-e*cos(E));
        E = E + delta_E;
    end
    
end

function Rx = Rx(alpha)

    Rx = [1 0 0 ;
            0 cos(alpha) -sin(alpha) ;
            0 sin(alpha) cos(alpha) ];
end

function Rz = Rz(alpha)

    Rz = [cos(alpha) -sin(alpha) 0 ; 
            sin(alpha) cos(alpha) 0 ;
             0 0 1 ];
end
