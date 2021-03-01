function moon = moon()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    solar_system.solar_system_parameters();
    moon = solar_system.body();
    moon.name = "moon";
    moon.mass = moon_mass;
    moon.diameter = moon_diameter;
    moon.a =moon_a;
    moon.i = moon_i;
    moon.e= moon_e;
    moon.an = moon_an;
    moon.pa = moon_pa;
    moon.l = moon_l;
end

