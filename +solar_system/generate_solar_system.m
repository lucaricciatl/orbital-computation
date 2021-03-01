function s = generate_solar_system()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
s = struct();

%---create bodies---%
s.sun = solar_system.sun();
s.earth = solar_system.earth();
s.moon = solar_system.moon();
s.mars = solar_system.mars();
s.jupiter = solar_system.jupiter();
% s.io = solar_system.io();
% s.callisto = solar_system.callisto();
% s.europa = solar_system.europa();
s.mercury = solar_system.mercury();
s.venus = solar_system.venus();
s.saturn = solar_system.saturn();
s.uranus = solar_system.uranus();
s.neptune = solar_system.neptune();
%---set attractors---%
s.earth.attractor = s.sun;
s.moon.attractor = s.earth;
s.mars.attractor = s.sun;
s.jupiter.attractor = s.sun;
% s.io.attractor = s.jupiter;
% s.callisto.attractor = jupiter;
% s.europa.attractor = jupiter;
s.mercury.attractor = s.sun;
s.venus.attractor = s.sun;
s.saturn.attractor = s.sun;
s.uranus.attractor = s.sun;
s.neptune.attractor = s.sun;

%initialize structure if necessary
% fn = fieldnames(s);
% num_struct = numel(fn);
% for k=1:num_struct
%     body = s.(fn{k});
%     body.compute_orbital_elements();  
% end

end

