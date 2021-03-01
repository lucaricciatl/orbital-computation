mtoau = 6.684589e-12;              %fattore di conversione metri unita astronomiche
sectodays = 0.0000115741;       %fattore di conversione secondi  giorni
G = 6.67408e-11;       %costante di gravitazione universale


%sun
sun.mass = 1.98855e30;

% earth
earth.a = 1.00000261;
earth.e = 0.01671123;
earth.I =0;%-0.00001531 * pi/180;
earth.L =100.46457166 * pi/180;
earth.om_bar =102.93768193 * pi/180;
earth.om_big =0;
earth.mass = 5.97219e24
earth.radius_of_influ = 9.246e8;
earth.radius = 6378388; %[m]
earth.mu = earth.mass*G; 

earth.a_dot =0.00000562;
earth.e_dot =-0.00004392;
earth.I_dot =0;%-0.01294668 * pi/180;
earth.L_dot =35999.37244981 * pi/180;
earth.om_bar_dot =0.032327364 * pi/180;
earth.om_big_dot =0;

% mars
mars.a = 1.52371034;
mars.e = 0.09339410;
mars.I =1.84969142 * pi/180;
mars.L=-4.55343205 * pi/180;
mars.om_bar=-23.94362959 * pi/180;
mars.om_big= 49.55953891 * pi/180;

mars.a_dot=0.00001847;
mars.e_dot=0.00007882;
mars.I_dot= -0.00813131 * pi/180;
mars.L_dot=19140.30268499 * pi/180;
mars.om_bar_dot=0.44441088 * pi/180;
mars.om_big_dot= -0.29257343 * pi/180;

%venus
venus.a = 0.72333566;
venus.e = 0.00677672;
venus.I = 0;%3.39467605 * pi/180;
venus.L= 181.97909950 * pi/180;
venus.om_bar=131.60246718 * pi/180;
venus.om_big= 76.67984255 * pi/180;

venus.a_dot= 0.00000390;
venus.e_dot=-0.00004107;
venus.I_dot=0;%-0.00078890 * pi/180;
venus.L_dot=58517.81538729 * pi/180;
venus.om_bar_dot=0.00268329 * pi/180;
venus.om_big_dot =-0.27769418 * pi/180;

%%mercury
% mercury.a = 0.38709927;
% mercury.e = 0.20563593;
% mercury.I = 7.00497902;
% mercury.L = 252.25032350;
% mercury.om_bar=77.45779628;
% mercury.om_big= 48.33076593;
% 
% mercury.a_dot
% mercury.e_dot
% mercury.I_dot
% mercury.L_dot
% mercury.om_bar_dot
% mercury.om_big_dot
