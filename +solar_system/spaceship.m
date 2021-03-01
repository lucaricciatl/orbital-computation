function spaceship = spaceship()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    solar_system.solar_system_parameters();
    spaceship = solar_system.body();
    spaceship.name = "earth";
    spaceship.e = earth_e;
    spaceship.a = earth_a;
    spaceship.i = earth_i;
    spaceship.an = earth_an;
    spaceship.pa = earth_pa;
    spaceship.mass = earth_mass;
    spaceship.diameter = earth_diameter;
    spaceship.l = earth_l;

    spaceship.a_dot = earth_a_dot;
    spaceship.l_dot = earth_l_dot;
    spaceship.e_dot = earth_e_dot ;
    spaceship.i_dot = earth_i_dot;
    spaceship.an_dot = earth_pa_dot;
    spaceship.pa_dot = earth_an_dot;
end

