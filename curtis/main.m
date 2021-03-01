%% Preparations
clc
clear all
close all
tic

%addpath(genpath("+curtis"));
%addpath(genpath("Dawn Computations"));

%% Constants
%Sun gravitational parameter
global mu
mu = 1.327565122000000e+11; %[km^3/s^2]

colors = ["g"          %green
          "m"          %magenta
          "b"          %blue
          "r"          %red
          "#A2142F"    %darker red
          "#7E2F8E"    %purple
          "#4DBEEE"    %darker cyan
          "c"          %(bright) cyan
          "#D95319"    %orange
          "#77AC30"    %darker green
          "#EDB120"    %ochre
          "#D95319"];  %orange, not visible due to Sun orbit dimensions
      
[xx,yy,zz] = sphere(5);
% Ordine distanza crescente dal sole:
%                1 = Sun
%                2 = Venere
%                3 = Earth
%                4 = Giove
%                5 = Io

%Masse 
masses = 10^24 * [1989100   %Sole      %[kg]
                  4.86732   %Venere
                  5.97219   %Terra      
                  1898.13   %Giove
                  0.089];   %Io

%Raggio medio               
radii = [695508      %Sole              %[km]
         6051.8      %Venere
         6371        %Terra
         69911       %Giove 
         1821.5];    %Io

%Distanza media in km 
% Anche se il valore del sole non torna (??)
% Abbiamo consideraTO Io alla stessa distanza dal sole di Giove 
distances = [413690250       %Sole   %[km]
             108209475       %Venere
             149598262       %Terra
             778340821       %Giove
             778340821];     %Io
    
G = 6.6742e-20; %[km^3/kg/s^2]  %Costante di gravitazione universale

%% Data

%Distances from the Sun [km]
Venere_to_Sun = distances(2);
Earth_to_Sun = distances(3); 
Io_to_Sun = distances(5);

%Masses [Kg]
Sun_mass = masses(1);
Venere_mass = masses(2);
Earth_mass = masses(3);
Giove_mass = masses(4);
Io_mass = masses(5);

%SOI: (m_planet/m_Sun)^(2/5) * distance_from_Sun RAGGIO SFERA DI INFLUENZA km]
Venere_SOI = (masses(2)/Sun_mass)^(2/5)*distances(2);
Earth_SOI = (masses(3)/Sun_mass)^(2/5)*distances(3); 
Giove_SOI = (masses(4)/Sun_mass)^(2/5)*distances(4); 
Io_SOI = (masses(5)/Sun_mass)^(2/5)*distances(5);

%Parking orbits (not considering body radius)  
Epark_radius = 200; %[km]
Epark_inclination = 0; %[rad]

% NON SAPPIAMO COSA SONO
%_hamo raggio dell'inzio del flyby in ingresso al pianeta
%_lamo raggio della fine del flyby in uscita al pianeta
% Venere_hamo = ; %[km]
% Venere_lamo = ; %[km]
% Terra_hamo = ; %[km]
% Terra_lamo = ; %[km]

%% Delta v vectors
%VETTORE CONTENENTE TUTTI I DELTAV DEL VIAGGIO
delta_v = zeros(10);
% delta_v(1) escape from Earth
% delta_v(2) capture on Venere
% delta_v(3) change orbit on Venere
% delta_v(4) escape from Venere   
% delta_v(5) capture on Terra   
% delta_v(6) change orbit on Terra                             
% delta_v(7) escape from Earth                             
% delta_v(8) capture on Giove
% delta_v(9) escape from Giove
% delta_v(10) trasferimento su Io 

%DELTAV PER CAMBIARE DA ORBITA INTORNO AL PIANETA A IPERBOLE DI FUGA (???)
delta_change = zeros(3);
% delta_change(1) change of plane on Earth ???
% delta_change(2) change of plane after Mars ????
% delta_change(3) change of plane on Vesta ???

%% Planets configurations
%{
    The mission is subdivided in four parts:
    - Earth to Venere;
    - Venere to Earth;
    - Earth to Jupiter;
    - Jupiter to Io.
    Variables for each part have the same name but different suffix.
    
% Fly-by considerato istantaneo
    Planet configurations have the following numerical indications:
    - 0: at Earth departure (to Venere) -> 07/09/2029
    - 1: at Venere arrival/departure (from Venere, to Earth) -> 15/09/2030
    - 2: at Earth arrival/departure (from Earth, to Jupiter) -> 15/04/2032
    - 3: at Jupiter arrival (from Earth) -> 15/10/2034
    - 4: at Jupiter departure (to io) -> ????
    - 5: at Io arrival (from Jupiter) -> ????
%}
%   coe  - vector of orbital elements [h e omega i w theta a]
%Earth_r0 POSIZIONE della terra AL TEMPO 0 (della tabella precedente)
%Earth_v0 VELOCITA' della terra AL TEMPO 0 (della tabella precedente)

global Earth_r0
global Earth_r1
global Earth_r2
global Venere_r0
global Venere_r1
global Venere_r2
global Giove_r2
global Giove_r3

%Earth
[Earth_coe0, Earth_r0, Earth_v0, ~] =...
                        planet_elements_and_sv(3,2029,9,7,0,0,0);
[Earth_coe1, Earth_r1, Earth_v1, ~] =...
                        planet_elements_and_sv(3,2030,9,15,0,0,0);
[Earth_coe2, Earth_r2, Earth_v2, ~] =...
                        planet_elements_and_sv(3,2032,4,15,0,0,0);

%Venere
[Venere_coe0, Venere_r0, Venere_v0, ~] =...
                        planet_elements_and_sv(2,2029,9,7,0,0,0);
[Venere_coe1, Venere_r1, Venere_v1, ~] =...
                        planet_elements_and_sv(2,2030,9,15,0,0,0);
[Venere_coe2, Venere_r2, Venere_v2, ~] =...
                        planet_elements_and_sv(2,2032,4,15,0,0,0);

%Giove
[Giove_coe2, Giove_r2, Giove_v2, ~] =...
                        planet_elements_and_sv(4,2032,4,15,0,0,0);
[Giove_coe3, Giove_r3, Giove_v3, ~] =...
                        planet_elements_and_sv(4,2034,10,15,0,0,0);
                    

% DOBBIAMO DECIDERE LA DATA DI PARTENZA
% [Giove_coe4, Giove_r4, Giove_v4, ~] =...
%                         planet_elements_and_sv(4,2009,2,17,0,0,0);
% [Giove_coe5, Giove_r5, Giove_v5, ~] =...
%                         planet_elements_and_sv(4,2009,2,17,0,0,0);

%Io
% DOBBIAMO DECIDERE LA DATA DI PARTENZA
% [Io_coe4, Io_r4, Io_v4, ~] =...
%                         planet_elements_and_sv(11,2012,9,5,0,0,0);  
% [Io_coe5, Io_r5, Io_v5, ~] =...
%                         planet_elements_and_sv(11,2015,3,5,0,0,0);  

%% Interplanetary plot
%% Earth - Venere travel (0 --> 1)
if exist('figure2') == 0  
    figure()
else
    figure2()
end

title("Interplanetary orbits")
hold on

%Interplanetary orbit
fprintf('\n\n EARTH TO VENERE \n\n')
%   [dep_r,dep_v,arr_r,arr_v,flight,orb_oe] = gen_orbit(...)
%   returns the starting and final position dep_r, arr_r
%   of the spacecraft (corresponding to the main bodies position
%   due to the use of the patched conics method), along with its
%   starting and final velocities dep_v, arr_v, the time of flight
%   (flight) and the orbital elements orb_oe of the trajectory.

[body_pos1, sp_v1, body_posf1, sp_vf1,tof1, orb_elem1] = ...
                gen_orbit(3,2,[2029 9 7 0 0 0],[2030 9 15 0 0 0],0);    
EM_orbit = intpl_orbit(tof1,Earth_r0,sp_v1); 
% body_pos1 POSIZIONE NAVICELLA AL TEMPO 1
%tof1 delta t per arrivare a destinazione  
%
                                                                      
%Planet orbits
plot_orbit(3,2029)
plot_orbit(2,2030)

%Planet positions
plot3(Earth_r0(1),Earth_r0(2),Earth_r0(3),'o','Color',colors(3))
plot3(Earth_r1(1),Earth_r1(2),Earth_r1(3),'x','Color',colors(3))
plot3(Venere_r0(1),Venere_r0(2),Venere_r0(3),'o','Color',colors(4))
plot3(Venere_r1(1),Venere_r1(2),Venere_r1(3),'x','Color',colors(4))

xlabel('x')
ylabel('y')
zlabel('z')
view(-100000,100000)
grid

%% Venere - Earth travel (1 --> 2) (FLY-BY)
%SIAMO SU VENERE!!!
%Interplanetary orbit
fprintf('\n\n VENERE to EARTH \n\n')
%   [dep_r,dep_v,arr_r,arr_v,flight,orb_oe] = gen_orbit(...)
%   returns the starting and final position dep_r, arr_r
%   of the spacecraft (corresponding to the main bodies position
%   due to the use of the patched conics method), along with its
%   starting and final velocities dep_v, arr_v, the time of flight
%   (flight) and the orbital elements orb_oe of the trajectory.

% [body_pos2, sp_v2, body_posf2, sp_vf2, tof2, orb_elem2] = ...
%                 gen_orbit(4,10,[2009 2 17 0 0 0],[2011 7 16 0 0 0],0);
%            
% [dep_r, dep_v, arr_r, arr_v, flight, orb_oe]
%                1 = Sun
%                2 = Venere
%                3 = Earth
%                4 = Giove
%                5 = Io
% QUESTA PARTE NON L'HO CAPITA!!! CREDO SERVA PER IL FLYBY
[body_pos21, sp_todawn, body_posf21, sp_vf21, tof21, orb_elem21] = ...
                gen_orbit(2,3,[2030 9 15 0 0 0],[2032 4 15 0 0 0],1);
% flag = 1 dal pianeta di partenza all'added point
% flag = 2 dall'added point al pianeta di arrivo
[body_pos22, sp_fromdawn, body_posf22, sp_vf22, tof22, orb_elem22] = ...
                gen_orbit(2,3,[2030 9 15 0 0 0],[2032 4 15 0 0 0],1);

%   arr_days  - days needed for the spacecraft to arrive at destination
%   r_init    - initial position of the spacecraft
%   v_init    - initial velocity of the spacecraft

MD_orbit = intpl_orbit(tof21,Earth_r1,sp_todawn);
DV_orbit = intpl_orbit(tof22,body_posf21,sp_fromdawn);
MV_orbit = [MD_orbit;DV_orbit];
            
% MV_orbit = intpl_orbit(tof2,Mars_r1,sp_v2);

%Planet orbits
plot_orbit(2,2030)
plot_orbit(3,2032)

%Planet positions
plot3(Venere_r1(1),Venere_r1(2),Venere_r1(3),'o','Color',colors(4))
plot3(Earth_r1(1),Earth_r1(2),Earth_r1(3),'o','Color',colors(3))
plot3(Earth_r2(1),Earth_r2(2),Earth_r2(3),'x','Color',colors(3))

% TENTATIVO DI GRAFICO UGUALE A QUELLO PRECEDENTE

% [body_pos1, sp_v1, body_posf1, sp_vf1,tof1, orb_elem1] = ...
%                 gen_orbit(2,3,[2030 9 15 0 0 0],[2032 4 15 0 0 0],0);    
% EM_orbit = intpl_orbit(tof1,Earth_r1,sp_v1);                            
% %tof1 delta t per arrivare a destinazione  
% 
% %
%                                                                       
% %Planet orbits
% plot_orbit(2,2020)
% plot_orbit(3,2032)
% 
% %Planet positions
% plot3(Earth_r1(1),Earth_r1(2),Earth_r1(3),'o','Color',colors(3))
% plot3(Earth_r2(1),Earth_r2(2),Earth_r2(3),'x','Color',colors(3))
% plot3(Venere_r1(1),Venere_r1(2),Venere_r1(3),'o','Color',colors(4))
% plot3(Venere_r2(1),Venere_r2(2),Venere_r2(3),'x','Color',colors(4))
% % %% Earth - Giove travel (2 --> 3)
% 
% %Interplanetary orbit
% fprintf('\n\n EARTH TO GIOVE \n\n')
% 
% %                1 = Sun
% %                2 = Venere
% %                3 = Earth
% %                4 = Giove
% %                5 = Io
% 
% % [body_pos3, sp_v3, body_posf3, sp_vf3, tof3, orb_elem3] = ...
% %                 gen_orbit(3,4,[2032 4 15 0 0 0],[2034 10 15 0 0 0],0);
% % VC_orbit = intpl_orbit(tof3,Vesta_r3,sp_v3);
% % 
% % %Planet orbits
% % plot_orbit(11,2015)
% % 
% % %Planet positions
% % plot3(Vesta_r3(1),Vesta_r3(2),Vesta_r3(3),'d','Color',colors(10))
% % plot3(Vesta_r4(1),Vesta_r4(2),Vesta_r4(3),'+','Color',colors(10))
% % plot3(Ceres_r3(1),Ceres_r3(2),Ceres_r4(3),'d','Color',colors(11))
% % plot3(Ceres_r4(1),Ceres_r4(2),Ceres_r4(3),'+','Color',colors(11))
% 
% [body_pos31, sp_todawn1, body_posf31, sp_vf31, tof31, orb_elem31] = ...
%                 gen_orbit(3,4,[2032 4 15 0 0 0],[2034 10 15 0 0 0],1);
% [body_pos32, sp_fromdawn1, body_posf32, sp_vf32, tof32, orb_elem32] = ...
%                 gen_orbit(3,4,[2032 4 15 0 0 0],[2034 10 15 0 0 0],2);
% % NON SAPPIAMO COSA SIA IL FLAG (ATTIVATO CON VALORI 1 E 2)
% %   arr_days  - days needed for the spacecraft to arrive at destination
% %   r_init    - initial position of the spacecraft
% %   v_init    - initial velocity of the spacecraft
% 
% MD_orbit = intpl_orbit(tof31,Giove_r2,sp_todawn1);
% DV_orbit = intpl_orbit(tof32,body_posf31,sp_fromdawn1);
% MV_orbit = [MD_orbit;DV_orbit];
%             
% % MV_orbit = intpl_orbit(tof2,Mars_r1,sp_v2);
% 
% %Planet orbits
% plot_orbit(3,2032)
% plot_orbit(4,2034)
% 
% %Planet positions
% plot3(Earth_r2(1),Earth_r2(2),Earth_r2(3),'o','Color',colors(3))
% plot3(Giove_r2(1),Giove_r2(2),Giove_r2(3),'o','Color',colors(5))
% plot3(Giove_r3(1),Giove_r3(2),Giove_r3(3),'x','Color',colors(5))