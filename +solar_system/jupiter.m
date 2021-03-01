function jupiter = jupiter()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    solar_system.solar_system_parameters();
    jupiter = solar_system.body();
    jupiter.name = "jupiter"; 
    jupiter.mass = jupiter_mass;
    jupiter.diameter = jupiter_diameter;
    jupiter.a = jupiter_a;
    jupiter.i = jupiter_i;
    jupiter.e= jupiter_e;
    jupiter.an = jupiter_an;
    jupiter.pa = jupiter_pa;
    jupiter.l = jupiter_l;
    
    jupiter.a_dot = jupiter_a_dot;
    jupiter.l_dot = jupiter_l_dot;
    jupiter.e_dot = jupiter_e_dot ;
    jupiter.i_dot = jupiter_i_dot;
    jupiter.an_dot = jupiter_pa_dot;
    jupiter.pa_dot = jupiter_an_dot;
end

