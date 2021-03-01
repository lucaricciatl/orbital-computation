function mercury = mercury()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    solar_system.solar_system_parameters();
    mercury = solar_system.body();
    mercury.name = "mercury";
    mercury.mass = mercury_mass;
    mercury.diameter = mercury_diameter;
    mercury.a =mercury_a;
    mercury.i = mercury_i;
    mercury.e= mercury_e;
    mercury.an = mercury_an;
    mercury.pa = mercury_pa;
    mercury.l = mercury_l;
    
    mercury.a_dot = mercury_a_dot;
    mercury.l_dot = mercury_l_dot;
    mercury.e_dot = mercury_e_dot ;
    mercury.i_dot = mercury_i_dot;
    mercury.an_dot = mercury_pa_dot;
    mercury.pa_dot = mercury_an_dot;
end

