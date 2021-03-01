function venus = venus()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    solar_system.solar_system_parameters();
    venus = solar_system.body();
    venus.mass = venus_mass;
    venus.name = "venus";
    venus.diameter = venus_diameter;
    venus.a = venus_a;
    venus.i = venus_i;
    venus.e= venus_e;
    venus.an = venus_an;
    venus.pa = venus_pa;
    venus.l = venus_l;
    
    venus.a_dot = venus_a_dot;
    venus.l_dot = venus_l_dot;
    venus.e_dot =venus_e_dot ;
    venus.i_dot = venus_i_dot;
    venus.an_dot = venus_pa_dot;
    venus.pa_dot = venus_an_dot;
end

