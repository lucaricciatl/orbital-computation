function saturn = saturn()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    solar_system.solar_system_parameters();
    saturn = solar_system.body();
    saturn.mass = saturn_mass;
    saturn.name = "saturn";
    saturn.diameter = saturn_diameter;
    saturn.a = saturn_a;
    saturn.i = saturn_i;
    saturn.e= saturn_e;
    saturn.an = saturn_a;
    saturn.pa = saturn_pa;
    saturn.l = saturn_l;
    
    saturn.a_dot = saturn_a_dot;
    saturn.l_dot = saturn_l_dot;
    saturn.e_dot = saturn_e_dot ;
    saturn.i_dot = saturn_i_dot;
    saturn.an_dot = saturn_pa_dot;
    saturn.pa_dot = saturn_an_dot;
end

