%% Charan Srikkanth - EE23B127 , Astha Chand - EE23B013
% This code contains our simulations of the endsem presentation for the course EE7500 on Digital Coding Metasurfaces



function far_field = get_far_field_for_code(code)
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
    figure;
    imagesc(code);
    colormap(gray);
    title('Metasurface Code');
    xlabel('Column Index');
    ylabel('Row Index');
end

%code = repmat([1 1 1 1 0 0 0 0;1 1 1 1 0 0 0 0;0 0 0 0 1 1 1 1;0 0 0 0 1 1 1 1;], 2, 1);
code = repmat([1 1 1 1 1 1 1 1;0 0 0 0 0 0 0 0;1 1 1 1 1 1 1 1;0 0 0 0 0 0 0 0;],2,1);
code = randi([0,1],config.M,config.N);
%code = repmat([1 0 1 0 1 0 1 0;0 1 0 1 0 1 0 1;1 0 1 0 1 0 1 0;0 1 0 1 0 1 0 1;],2,1);
%code = zeros(config.M,config.N);
output_field =get_far_field_for_code(code);

E_mag = abs(output_field);

% 2. Normalize and convert to dB for the color map
E_norm = E_mag / max(max(E_mag));  %Max does only on one axis we want for both
R_dB = mag2db(E_norm);

% Apply a floor of -40 dB so the color scale isn't skewed by infinite nulls
R_dB(R_dB < -40) = -40;

[Phi,Theta] = meshgrid(config.phi,config.theta);
%R = pow2db(abs(output_field));
X = E_mag.*sin(Theta).*cos(Phi);
Y = E_mag.*sin(Theta).*sin(Phi);
Z = E_mag.*cos(Theta);

fig = figure('Visible','off');
surf(X,Y,Z,R_dB);
shading interp;   % Smooths out the grid lines for a cleaner look
colormap jet;     % Standard heatmap colors
colorbar;         % Adds a reference scale on the side
% 6. Formatting the plot
axis equal;       % CRITICAL: Ensures the 3D shape isn't artificially stretched
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Far-Field Radiation Pattern');
view(45, 30);     % Adjusts the default camera angle


saveas(fig,'far_field_pattern_3d.png'); % Save the figure as a PNG file

fig2=figure('visible','off');
phi = linspace(0, 360, size(E_mag, 2));
theta = linspace(0, 90, size(E_mag, 1));
pcolor(phi, theta, E_mag);
shading interp;
xlabel('\phi (degrees)');
ylabel('\theta (degrees)');
colorbar;
title('E_{mag} vs \phi and \theta');
saveas(fig2,'far_field_pattern_2d.png'); % Save the figure as a PNG file

plot_code(code);
%% Verifying the convolution property of the array factor for 2 bit array
function far_field = get_far_field_for_2d_code(code)
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
    title(figname);
    view(45, 30);     % Adjusts the default camera angle
    filename = sprintf('%s.png', figname);
    saveas(fig,filename); % Save the figure as a PNG file
end

function code_plot = plot_2bit_code(code,figname)
    fig = figure('Visible','off');
    imagesc(code);
    colormap(gray(4)); % Use a 4-level grayscale colormap (white to black)
    colorbar('Ticks',0:3,'TickLabels',{'0','1','2','3'});
    caxis([0 3]);
    title(figname);
    xlabel('Column Index');
    ylabel('Row Index');
    filename = sprintf('%s.png', figname);
    saveas(fig,filename); % Save the figure as a PNG file
end

plus_code = zeros(config.M,config.N);
plus_code(4,:) = 2;
plus_code(5,:) = 2;
plus_code(:,4) = 2;
plus_code(:,5) = 2;
plot_2bit_code(plus_code,'Plus Code for 2-bit Array');
%code = repmat([1 0 1 0 1 0 1 0;0 1 0 1 0 1 0 1;1 0 1 0 1 0 1 0;0 1 0 1 0 1 0 1;],2,1);
far_field_2bit = get_far_field_for_2d_code(plus_code);
plot_far_field_2d(far_field_2bit,'Far-Field Pattern for 2-bit Array');

gradient_code = repmat(0:3,8,2);
plot_2bit_code(gradient_code,'Gradient Code for 2-bit Array');
far_field_2bit_gradient = get_far_field_for_2d_code(gradient_code);
plot_far_field_2d(far_field_2bit_gradient,'Far-Field Pattern for Gradient Code');

net_code = plus_code + gradient_code;
plot_2bit_code(net_code,'Net Code for 2-bit Array');
far_field_2bit_net = get_far_field_for_2d_code(net_code);
plot_far_field_2d(far_field_2bit_net,'Far-Field Pattern for Net Code');
