function io = io()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    solar_system.solar_system_parameters();
    io = solar_system.body();
    io.name = "io";
    io.mass = io_mass;
    io.diameter = io_diameter;
    io.a = io_a;
    io.i = io_i;
    io.e= io_e;
    io.an = io_an;
    io.pa = io_pa;
    io.l = io_l;
end

