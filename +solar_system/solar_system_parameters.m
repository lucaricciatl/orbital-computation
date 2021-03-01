
%-----------------------------------------------------------------------%
%                                                                       %
%                       SOLAR SYSTEM ORBITAL PARAM                      %
%                                                                       %
%-----------------------------------------------------------------------%
% Sono richiesti sei parametri per specificare l'orbita kepleriana di un
% corpo. Ad esempio, i tre numeri che descrivono la posizione iniziale del
% corpo e i tre valori per la sua velocità descriveranno un'orbita unica
% che può essere calcolata sia in avanti che indietro. Tuttavia, i
% parametri solitamente utilizzati sono leggermente diversi.

%---
%p0 = x0,y0,z0
%v0 = vx0,vy0,vz0
%---

%-----------------------------------------------------------------------%
%ORBITAL PARAMETERS DESCRIPTION
%-----------------------------------------------------------------------%
%here the animation : https://www.youtube.com/watch?v=QZrYaKwZwhI
%a -> major semiaxis (Km)
%e -> eccentricity 
%i -> inclination [deg]
%an -> ascending node longitude [deg] (big omega)
%pa -> Argument of periapsis [deg] (little omega)
%L -> mean longitude (L)


%-----------------------------------------------------------------------%
%                               CONSTANTS
%-----------------------------------------------------------------------%
G = 6.67*10^-11; %[m^3/(kg*s^2)]
AU = 1/149597870700; %[m/au]
DR = 2*pi/360; %%deg to radiants conversion units
%-----------------------------------------------------------------------%
%                               SUN
%-----------------------------------------------------------------------%
%used as main attractor
sun_mass = 1.989 * 10^30;
sun_diameter = 1392700 * 10^3 * AU; %%scaled by 20 for visualization
sun_x0 = 0;
sun_y0 = 0;
sun_z0 = 0;
sun_vx0 = 0;
sun_vy0 = 0;
sun_vz0 = 0;
sun_i = 0;
sun_a = 0;
sun_e = 0;
sun_an = 0 * DR;
sun_pa =0 * DR;
sun_l = 0; %%%MISSING DATAS

%-----------------------------------------------------------------------%
%                           EARTH SYSTEM
%-----------------------------------------------------------------------%
%---EARTH---
earth_a = 149597887*10^3 * AU;
earth_e = 0.01671022;
earth_i = -0.0000153 * DR;
earth_an = 0 * DR;
earth_pa = 102.20783 * DR;
earth_l = 100.94719*DR ; %100.46435; 
earth_mass = 5.9726*10^24;
earth_diameter = 12745594 * AU;
earth_a_dot =0.00000562 ;
earth_e_dot =-0.00004392 ;
earth_i_dot =-0.01294668  * DR;
earth_l_dot = 35999.37244981  * DR;
earth_pa_dot =0.032327364  * DR;
earth_an_dot =0 * DR;
earth_influence_radius = 1500000*AU;
%--MOON--
moon_a = 384400*10^3 * AU;
moon_e = 0.0549 * DR;
moon_i = 5.145 * DR;
moon_an = 0;% MISSING DATAS;
moon_pa = 0;%MISSING DATAS;
moon_l = 0;%MISSING DATAS;
moon_mass = 7.342*10^22;
moon_diameter = 3.476 * 10^3 * AU;
%-----------------------------------------------------------------------%
%                               MARS SYSTEM
%-----------------------------------------------------------------------%
%-mars
mars_a = 227939200 * 10^3 * AU;
mars_e = 0.09341233;
mars_i = 1.85 * DR;
mars_an = 49.57854*DR;
mars_pa = -23.943*DR ; 
mars_l = -4.553 * DR;
mars_mass = 6.4185*10^23;
mars_diameter = 6804.9 * 10^3 * AU;
mars_a_dot = 0.00001847 ;
mars_e_dot =- 0.00007288 ;
mars_i_dot = -0.008131  * DR;
mars_l_dot = 19140.3  * DR;
mars_pa_dot = 0.4444  * DR;
mars_an_dot = -0.29257343 * DR;
%-----------------------------------------------------------------------%
%                               JUPITER SYSTEM
%-----------------------------------------------------------------------%
%---JUPITER---
jupiter_a = 778412027*10^3 * AU;
jupiter_e = 0.04839266;
jupiter_i = 1.31 * DR;
jupiter_an = 100.4739 * DR;
jupiter_pa = 14.728 * DR;
jupiter_l = 34.40438 * DR ;
jupiter_mass = 1.89819*10^27;
jupiter_diameter = 142984 * 10^3 * AU;
jupiter_a_dot = 0.00011607 ;
jupiter_e_dot =- 0.00013253 ;
jupiter_i_dot = -0.008137  * DR;
jupiter_l_dot = 3034.74612  * DR;
jupiter_pa_dot = 0.21252  * DR;
jupiter_an_dot = -0.20469 * DR;
%--IO--
io_a = 421700*10^3 * AU;
io_e = 0.0041;
io_i = 0.05 * DR ;
io_an = 0;%MISSING DATAS
io_pa = 0;%MISSING DATAS
io_l = 0;%MISSING DATAS
io_mass =  8.931938*10^22;
io_diameter = 2*1821.6*10^3 * AU;

%--EUROPA--
europa_a = 671034*10^3 * AU;
europa_e = 0.0094;
europa_i = 0.47 * DR ;
europa_an =0; %MISSING DATAS
europa_pa =0 ;%MISSING DATAS
europa_l =0 ;%MISSING DATAS
europa_mass = 4.80*10^22;
europa_diameter = 3121.6 * 10^3 * AU;

%--CALLISTO--
% callisto_x0 = 
% callisto_y0 =
% callisto_z0 =
% callisto_vx0 =
% callisto_vy0 =
% callisto_vz0 =
callisto_a = 1882700*10^3 * AU;
callisto_e = 0.0074 ;
callisto_i = 2.21*10^3 * DR ;
callisto_an = 0 ;% MISSING DATAS;
callisto_pa =0; %MISSING DATAS;
callisto_l =0; %MISSING DATAS;
callisto_mass = 1.0759*10^23;
callisto_diameter = 4820.6 * 10^6 * AU; 
%-----------------------------------------------------------------------%
%                           MERCURY SYSTEM
%-----------------------------------------------------------------------%
%---MERCURY---
mercury_a = 5.791*10^10 * AU;
mercury_e= 0.2056 ;
mercury_i = 7.01 * DR ;
mercury_an = 48.33167 * DR ;
mercury_pa = 77.45645 *DR ;
mercury_l = 252.25084*DR;
mercury_mass = 3*3011*10^23;
mercury_diameter = 4879.4*10^3 * AU;
mercury_a_dot = 0.00000037 ;
mercury_e_dot = -0.0001906 ;
mercury_i_dot = -0.005947  * DR;
mercury_l_dot = 149472.67  * DR;
mercury_pa_dot = 0.16047  * DR;
mercury_an_dot = -0.1253 * DR;
%-----------------------------------------------------------------------%
%                           VENUS SYSTEM
%-----------------------------------------------------------------------%
%-venus
venus_a = 1.0821*10^11 * AU;
venus_e = 0.0067;
venus_i =	3.39 * DR;
venus_an = 76.68069 * DR;
venus_pa = 131.53298 *DR ;%54.8522* DR;
venus_l = 181.97973*DR;
venus_mass = 4.8675*10^24;
venus_diameter = 12103.6 *10^6 * AU;
venus_a_dot = 0.0000039 ;
venus_e_dot =- 0.00004107 ;
venus_i_dot = -0.00078890  * DR;
venus_l_dot =58517.815  * DR;
venus_pa_dot = 0.002683  * DR;
venus_an_dot = -0.2777 * DR;
%-----------------------------------------------------------------------%
%                               SATURN SYSTEM
%-----------------------------------------------------------------------%
%--SATURN---
saturn_a = 1433530000*10^3 * AU;
saturn_e = 0.0565;
saturn_i = 2.485* DR;
saturn_an = 113.66 * DR;
saturn_pa = 92.598* DR;
saturn_l = 49.95 *DR ;% 49.94432* DR; 
saturn_mass = 5.6834*10^26;
saturn_diameter = 120536*10^3 * AU;
saturn_a_dot = -0.0012506 ;
saturn_e_dot = -0.00050991 ;
saturn_i_dot = -0.001936  * DR;
saturn_l_dot = 1222.9  * DR;
saturn_pa_dot = -0.4187  * DR;
saturn_an_dot = -0.28867794 * DR;
%-----------------------------------------------------------------------%
%                           URANUS SYSTEM
%-----------------------------------------------------------------------%
%---URANUS---
uranus_a = 2872.46*10^9 * AU;
uranus_e = 0.0457;
uranus_i = 0.76986* DR;
uranus_an = 74.22988* DR;
uranus_pa = 170.96424 *DR ;% 96.541* DR;
uranus_l = 313.23218* DR;
uranus_mass = 86813 * 10^24;
uranus_diameter = 51118 * 10^3 * AU;
uranus_a_dot = -0.00196176 ;
uranus_e_dot = -0.000043971 ;
uranus_i_dot = -0.0024293  * DR;
uranus_l_dot = 428.48  * DR;
uranus_pa_dot = -0.408 * DR;
uranus_an_dot = -0.042 * DR;
%-----------------------------------------------------------------------%
%                           NEPTUNE SYSTEM
%-----------------------------------------------------------------------%
%---NEPTUNE---
neptune_a = 4498252900*10^3 * AU;
neptune_e = 0.00858587;
neptune_i = 1.77* DR;
neptune_an = 131.72169* DR;
neptune_pa = 44.97135 *DR ; %273.24966* DR;
neptune_l = 304.88003*DR;
neptune_mass = 1.0243*10^26 ;
neptune_diameter = 48681*10^3 * AU;
neptune_a_dot = -0.00026291 ;
neptune_e_dot = -0.000051051 ;
neptune_i_dot = -0.00035372  * DR;
neptune_l_dot = 281.45  * DR;
neptune_pa_dot = -0.322 * DR;
neptune_an_dot = -0.011834 * DR;
