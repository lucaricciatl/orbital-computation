function uranus = uranus()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    solar_system.solar_system_parameters();
    uranus = solar_system.body();
    uranus.mass = uranus_mass;
    uranus.name = "uranus";
    uranus.diameter = uranus_diameter;
    uranus.a = uranus_a;
    uranus.i = uranus_i;
    uranus.e= uranus_e;
    uranus.an = uranus_an;
    uranus.pa = uranus_pa;
    uranus.l = uranus_l;
    
    uranus.a_dot = uranus_a_dot;
    uranus.l_dot = uranus_l_dot;
    uranus.e_dot = uranus_e_dot ;
    uranus.i_dot = saturn_i_dot;
    uranus.an_dot = uranus_pa_dot;
    uranus.pa_dot = uranus_an_dot;
end

