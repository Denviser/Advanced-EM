classdef config
    % This class stores global constants for the array simulation
    properties (Constant)
        fc = 8.6e9;
        M = 8;
        N = 8;
        c= 3e8;
        k_0 = 2*pi*config.fc/config.c;
        lambda = config.c/config.fc;
        Dx = config.lambda;  %Set it to lambda 
        Dy = config.lambda;  %Set it to lambda
        n_points_theta = 1000;
        n_points_phi = 1000;
        theta = linspace(0,pi/2,config.n_points_theta);
        phi = linspace(0,2*pi,config.n_points_phi);
    end
end