function neptune = neptune()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    solar_system.solar_system_parameters();
    neptune = solar_system.body();
    neptune.name = "neptune";
    neptune.mass = neptune_mass;
    neptune.diameter = neptune_diameter;
    neptune.a = neptune_a;
    neptune.i = neptune_i;
    neptune.e= neptune_e;
    neptune.an = neptune_an;
    neptune.pa = neptune_pa;
    neptune.l = neptune_l;
    
    neptune.a_dot = neptune_a_dot;
    neptune.l_dot = neptune_l_dot;
    neptune.e_dot = neptune_e_dot ;
    neptune.i_dot = neptune_i_dot;
    neptune.an_dot = neptune_pa_dot;
    neptune.pa_dot = neptune_an_dot;
end

