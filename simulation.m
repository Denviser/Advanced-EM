%% Charan Srikkanth - EE23B127 , Astha Chand - EE23B013
% The global constants used in this file are defined in the config.m file
% This code contains our simulations of the endsem presentation for the course EE7500 on Digital Coding Metasurfaces.
% It contains two sections the first section you can just create an 2D array of 0 and 1 (1bit metasurface code) and the plots for the code,3D far field pattern and 2D far field pattern are saved in ./plots/one_bit.
% The second section is to show the convolution property of a 2 bit code(explanation can be found in report). Similar to the one bit it takes in the code and plots the 3D far field pattern and 2D far field pattern. The plots are saved in ./plots/two_bit
% Make sure to to have the corresponding folders ./plots/one_bit and ./plots/two_bit created before running the code.

function far_field = get_far_field_for_code(code)
    % This function takes in the code for the surface(0 or 1) and returns a matrix of far field values for each theta and phi
    far_field = zeros(config.n_points_theta,config.n_points_phi);
    for particular_theta_index=1:config.n_points_theta
        for i=1:config.M
            for j=1:config.N
                if code(i,j) == 0
                    code_phase = 0;
                else
                    code_phase = pi;
                end
                
                element_factor = cos(config.theta(particular_theta_index)); % Cosine element factor for a simple dipole
                far_field(particular_theta_index,:) = far_field(particular_theta_index,:) + element_factor.*exp(1j*code_phase).*exp(1j*config.k_0.*sin(config.theta(particular_theta_index)).*(config.Dx*(i-0.5).*cos(config.phi)+config.Dy*(j-0.5).*sin(config.phi)));
            
            end
        end
    end
end


function code_plot = plot_code(code)
    %This function just takes in the metasurface code and plots it.
    fig = figure('Visible','off');
    imagesc(code);
    colormap(gray);
    title('Metasurface Code');
    xlabel('Column Index');
    ylabel('Row Index');
    c=colorbar;
    clim([0 1]);
    c.Label.String = 'Code Value';
    saveas(fig,'./plots/one_bit/code.png');
end

%Uncomment the code according to which code you want to use

code = repmat([1 1 0 0 1 1 0 0;1 1 0 0 1 1 0 0;0 0 1 1 0 0 1 1;0 0 1 1 0 0 1 1;], 2, 1);
%code = repmat([1 1 1 1 1 1 1 1;0 0 0 0 0 0 0 0;1 1 1 1 1 1 1 1;0 0 0 0 0 0 0 0;],2,1);
%code = randi([0,1],config.M,config.N);
%code = repmat([1 0 1 0 1 0 1 0;0 1 0 1 0 1 0 1;1 0 1 0 1 0 1 0;0 1 0 1 0 1 0 1;],2,1);
%code = ones(config.M,config.N);


output_field =get_far_field_for_code(code);

E_mag = abs(output_field);

% 2. Normalize and convert to dB for the color map
E_norm = E_mag / max(max(E_mag));  %Max does only on one axis we want for both
R_dB = mag2db(E_norm); % The reason we create this R_dB is to colour the graph nicely according to it.

% Apply a floor of -40 dB so the color scale isn't skewed by infinite nulls
R_dB(R_dB < -40) = -40;

[Phi,Theta] = meshgrid(config.phi,config.theta);
X = E_mag.*sin(Theta).*cos(Phi);
Y = E_mag.*sin(Theta).*sin(Phi);
Z = E_mag.*cos(Theta);

fig = figure('Visible','off');
surf(X,Y,Z,R_dB);
shading interp;   % Smooths out the grid lines for a cleaner look
colormap jet;     % Standard heatmap colors
c=colorbar;         % Adds a reference scale on the side
c.Label.String = 'Magnitude(dB)';
c.Label.FontSize = 12;
c.Label.FontWeight = 'bold';
c.Label.Color = 'black';
% 6. Formatting the plot
axis equal;       % CRITICAL: Ensures the 3D shape isn't artificially stretched
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Far-Field Radiation Pattern');
view(45, 30);     % Adjusts the default camera angle


saveas(fig,'./plots/one_bit/far_field_pattern_3d.png'); % Save the figure as a PNG file

plot_far_field_2d(output_field,'one_bit/far_field_pattern_2d');

plot_code(code);
%% Verifying the convolution property of the array factor for 2 bit array
function far_field = get_far_field_for_2bit_code(code)
    % This function takes in the 2_bit code for the surface(0,1,2,3) and returns a matrix of far field values for each theta and phi
    far_field = zeros(config.n_points_theta,config.n_points_phi);
    for particular_theta_index=1:config.n_points_theta
        for i=1:config.M
            for j=1:config.N
                if code(i,j) == 0
                    code_phase = 0;
                elseif code(i,j) == 1
                    code_phase = pi/2;
                elseif code(i,j) == 2
                    code_phase = pi;
                else
                    code_phase = 3*pi/2;
                end
                
                element_factor = cos(config.theta(particular_theta_index)); % Cosine element factor for a simple dipole
                far_field(particular_theta_index,:) = far_field(particular_theta_index,:) + element_factor.*exp(1j*code_phase).*exp(1j*config.k_0.*sin(config.theta(particular_theta_index)).*(config.Dx*(i-0.5).*cos(config.phi)+config.Dy*(j-0.5).*sin(config.phi)));
            
            end
        end
    end
end

function plot_far_field_2d(far_field,figname)
    % This function takes in the far_field array generate from get_far_field function and plots the 2D far field pattern
    E_mag = abs(far_field);
    %E_norm = E_mag / max(max(E_mag));  % Normalize for color mapping
    R_dB = mag2db(E_mag);
    R_dB(R_dB < -40) = -40;  % Floor at -40 dB

    fig2=figure('visible','off');
    phi = linspace(0, 360, size(E_mag, 2));
    theta = linspace(0, 90, size(E_mag, 1));
    pcolor(phi, theta, R_dB);
    shading interp;
    xlabel('\phi (degrees)');
    ylabel('\theta (degrees)');
    c=colorbar;
    clim([-40 30]);
    c.Label.String = 'Magnitude(dB)';
    c.Label.FontSize = 12;
    c.Label.FontWeight = 'bold';
    c.Label.Color = 'black';
    title('E_{mag} vs \phi and \theta');
    saveas(fig2,['./plots/' figname '_2d.png']); % Save the figure as a PNG file
end

function plot_far_field_3d(far_field,figname)
    % This function takes in the far_field array generate from get_far_field function and plots the 2D far field pattern
    E_mag = abs(far_field);
    E_norm = E_mag / max(max(E_mag));  % Normalize for color mapping
    R_dB = mag2db(E_norm);
    R_dB(R_dB < -40) = -40;  % Floor at -40 dB

    [Phi,Theta] = meshgrid(config.phi,config.theta);
    X = E_mag.*sin(Theta).*cos(Phi);
    Y = E_mag.*sin(Theta).*sin(Phi);
    Z = E_mag.*cos(Theta);

    fig =figure('Visible','off');
    surf(X,Y,Z,R_dB);
    shading interp;   % Smooths out the grid lines for a cleaner look
    colormap jet;     % Standard heatmap colors
    colorbar;         % Adds a reference scale on the side
    % 6. Formatting the plot
    axis equal;       % CRITICAL: Ensures the 3D shape isn't artificially stretched
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title("3D Far-Field Radiation Pattern");
    view(45, 30);     % Adjusts the default camera angle
    saveas(fig,['./plots/' figname '.png']); % Save the figure as a PNG file
end

function code_plot = plot_2bit_code(code,figname)
    %This function just takes in the 2bit metasurface code and plots it.
    fig = figure('Visible','off');
    imagesc(code);
    colormap(gray(4)); % Use a 4-level grayscale colormap (white to black)
    colorbar('Ticks',0:3,'TickLabels',{'0','1','2','3'});
    caxis([0 3]);
    title("Two Bit Metasurface Code");
    xlabel('Column Index');
    ylabel('Row Index');
    saveas(fig,['./plots/' figname '.png']); % Save the figure as a PNG file
end

% Generating the plus code
plus_code = zeros(config.M,config.N);
plus_code(4,:) = 2;
plus_code(5,:) = 2;
plus_code(:,4) = 2;
plus_code(:,5) = 2;

plot_2bit_code(plus_code,'/two_bit/Plus-Code');
far_field_2bit = get_far_field_for_2bit_code(plus_code);
plot_far_field_3d(far_field_2bit,'/two_bit/Far-Field-Plus-Code');

% Generating the gradient code
gradient_code = repmat(0:3,8,2);
plot_2bit_code(gradient_code,'/two_bit/Gradient-Code');
far_field_2bit_gradient = get_far_field_for_2bit_code(gradient_code);
plot_far_field_3d(far_field_2bit_gradient,'/two_bit/Far-Field-Gradient-Code');

% Generating the net code
net_code = plus_code + gradient_code;
plot_2bit_code(net_code,'/two_bit/Net Code');
far_field_2bit_net = get_far_field_for_2bit_code(net_code);
plot_far_field_3d(far_field_2bit_net,'/two_bit/Far-Field-Net-Code');
