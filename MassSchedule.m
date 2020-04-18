function [u] = MassSchedule(t, g1, g2, g3, g4)
    % code to calculate the accelerations of the arms for the planar trrex
    
    u = NaN(4, 1);
    u(1) = g1(t);
    u(2) = g2(t);
    u(3) = g3(t);
    u(4) = g4(t);
    
end