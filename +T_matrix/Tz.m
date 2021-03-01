function Tz = Tz(alpha)
%TX Summary of this function goes here
%   Detailed explanation goes here
    Tz = [cos(alpha) -sin(alpha) 0 0; 
          sin(alpha) cos(alpha) 0 0;
          0 0 1 0;
          0 0 0 1];
end

