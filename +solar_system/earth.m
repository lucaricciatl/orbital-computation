function earth = earth()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    solar_system.solar_system_parameters();
    earth = solar_system.body();
    earth.name = "earth";
    earth.e = earth_e;
    earth.a = earth_a;
    earth.i = earth_i;
    earth.an = earth_an;
    earth.pa = earth_pa;
    earth.mass = earth_mass;
    earth.diameter = earth_diameter;
    earth.l = earth_l;
    earth.influence_radius = earth_influence_radius
    earth.a_dot = earth_a_dot;
    earth.l_dot = earth_l_dot;
    earth.e_dot = earth_e_dot ;
    earth.i_dot = earth_i_dot;
    earth.an_dot = earth_pa_dot;
    earth.pa_dot = earth_an_dot;
end

