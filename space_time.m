%% EE7500: Space-Time Coding Metasurface Simulation 
% This file contains our simulations for space time coding. We basically generate a 3D array for a space time code then plot 
% the phase,amplitude code for each harmonic. We also plot the field for each harmonic.
% Make sure to have the folders ./plots/space_time_amp,./plots/space_time_phase,./plots/space_time_2D_far_field created before running the code.
% You can change the code in the Main Simulation Parameters and Initialisation section to generate a 2 bit code.

%% 1. Helper Functions
function a_m = get_fourier_coeff(m, L, coding_sequence)
    % This function basically just takes in the time-coding sequence and returns the Fourier Coefficient of the m-th harmonic(which is the effective static code for that harmonic)
    % Here L is number of time slots and m is the index of the harmonic
    a_m = 0;
    for n = 1:L
        Gamma_n = coding_sequence(n); % Complex reflection Gamma^n_pq
        
        % MATLAB's sinc(x) is sin(pi*x)/(pi*x)
        sinc_term = sinc(m/L); 
        phase_term = exp(-1j * pi * m * (2*n - 1) / L);
        a_m = a_m + (Gamma_n / L) * sinc_term * phase_term;
    end
end

function [amp_harmonic_code,phase_harmonic_code] = get_harmonic_code(m,L,st_complex_code)
    % L is the number of time slots and m is the index of the harmonic 
    % This function returns the amplitude and phase code for the m-th harmonic
    amp_harmonic_code = zeros(config.M,config.N);   
    phase_harmonic_code = zeros(config.M,config.N);
    for p = 1:config.M
        for q = 1:config.N
           a_m_pq = get_fourier_coeff(m,L,squeeze(st_complex_code(p,q,:)));
           amp_harmonic_code(p,q) = abs(a_m_pq);
           phase_harmonic_code(p,q) = rad2deg(angle(a_m_pq)); %The phase is in degrees 
        end
    end
end

function F_m = get_harmonic_field(m, L, st_complex_code, config)
    % This function takes in the code for a particular harmonic and returns the field for that harmonic for all theta phi.
    F_m = zeros(config.n_points_theta, config.n_points_phi);
    [Phi, Theta] = meshgrid(config.phi, config.theta);
    
    for p = 1:config.M
        for q = 1:config.N
            % Extract time-coding sequence for element (p,q)
            seq = squeeze(st_complex_code(p, q, :));
            
            % Calculate Fourier Coefficient a^m
            a_m_pq = get_fourier_coeff(m, L, seq);
            
            % Element Factor: cos(theta)
            element_factor = cos(config.theta).'; 
            % Spatial Phase Term
            spatial_phase = exp(1j * config.k_0 * sin(Theta) .* ...
                (config.Dx*(p-1)*cos(Phi) + config.Dy*(q-1)*sin(Phi)));
            
            F_m = F_m + (element_factor .* a_m_pq .* spatial_phase);
        end
    end
end
%% 2. Main Simulation Parameters & Initialization
L = 8; 
harmonics = -3:3;
st_code = zeros(config.M, config.N, L);
st_complex = zeros(config.M, config.N, L);

% This generates the example space time code which is present in our report
for n = 1:L
    for p = 1:config.M
        for q = 1:config.N
            % Diagonal logic: state is 1 if (p == n) and (q == n), else 0
            % To match the 'moving diagonal' in image (b):
            if p == (L - n + 1)
                code = 1;
            else
                code = 0;
            end
            
            st_code(p, q, n) = code;
            % Map 1 to phase pi, 0 to phase 0 (1-bit coding)
            st_complex(p, q, n) = exp(1j * code * pi);
        end
    end
end

%Uncomment the below two lines to generate the fine beam control code.

% st_code = generate_special_st_code();
% st_complex = exp(1j .* pi/2 .* st_code);


%% 3. Calculate and Plot Harmonics
for m = harmonics
    fprintf('Processing Harmonic m = %d...\n', m);
    F_m = get_harmonic_field(m, L, st_complex, config);
    [amp_harmonic_code,phase_harmonic_code] = get_harmonic_code(m, L, st_complex);
    plot_amp_code(amp_harmonic_code, sprintf('m = %d', m));
    plot_phase_code(phase_harmonic_code, sprintf('m = %d', m));
    plot_far_field_2d(F_m, sprintf('m = %d', m),m);
end
%% 4. Generate Video of 1-bit Coding States
video_filename = './plots/spacetime_metasurface_1bit.mp4';
v = VideoWriter(video_filename, 'MPEG-4');
v.FrameRate = 4; 
open(v);
fig_vid = figure('Visible', 'off'); 
for n = 1:L
    imagesc(st_code(:,:,n));
    colormap(gray(2)); 
    colorbar('Ticks', 0:1, 'TickLabels', {'0','1'});
    clim([0 1]); 
    
    title(sprintf('1-bit Time Slot n = %d / %d', n, L));
    xlabel('Column (q)');
    ylabel('Row (p)');
    axis square;
    
    frame = getframe(fig_vid);
    writeVideo(v, frame);
    clf(fig_vid);
end
close(v);
fprintf('Video successfully saved: %s\n', video_filename);
%% Supporting Plot Function
function plot_far_field_3d(far_field, figname, config)
    % This function takes in the far_field array generate from get_far_field function and plots the 3D far field pattern
    E_mag = abs(far_field);
    E_norm = E_mag / max(max(E_mag));
    R_dB = mag2db(E_norm);
    R_dB(R_dB < -40) = -40;
    
    [Phi, Theta] = meshgrid(config.phi, config.theta);
    X = E_mag .* sin(Theta) .* cos(Phi);
    Y = E_mag .* sin(Theta) .* sin(Phi);
    Z = E_mag .* cos(Theta);
    
    fig = figure('Visible', 'off');
    surf(X, Y, Z, R_dB);
    shading interp;
    colormap jet;
    colorbar;
    axis equal;
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(figname);
    view(45, 30);
    saveas(fig, ['./plots/' figname '.png']);
end

function code_plot = plot_code(code,figname)
    % This function just plots the code
    fig = figure('Visible','off');
    imagesc(code);
    colormap(gray);
    title(figname);
    xlabel('Column Index');
    ylabel('Row Index');
    c=colorbar;
    c.Label.String = 'Code Value';
    saveas(fig,['./plots/' figname '.png']);
end

function output = plot_amp_code(code,figname)
    % This function plots the amplitude code for a particular harmonic
    fig = figure('Visible','off');
    imagesc(code);
    colormap(hot);
    clim([0 1]);
    title(['Amplitude harmonic code ' figname]);
    xlabel('Column Index');
    ylabel('Row Index');
    c=colorbar;
    c.Label.String = 'Amplitude';
    saveas(fig,['./plots/space_time_amp/' figname '.png']);
end

function output = plot_phase_code(phase_code,figname)
    % This function plots the phase code for a particular harmonic
    fig = figure('Visible','off');
    imagesc(phase_code);
    colormap(winter);
    clim([-180 180]);
    title(['Phase harmonic code ' figname]);
    xlabel('Column Index');
    ylabel('Row Index');
    c=colorbar;
    c.Label.String = 'Code Value';
    saveas(fig,['./plots/space_time_phase/' figname '.png']);
end

function plot_far_field_2d(far_field,figname,m)
    % This function takes in the far_field array generate from get_far_field function and plots the 2D far field pattern for a particular harmonic
    E_mag = abs(far_field);
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
    colormap(parula)
    clim([-40 30]);
    c.Label.String = 'Magnitude(dB)';
    c.Label.FontSize = 12;
    c.Label.FontWeight = 'bold';
    c.Label.Color = 'black';
    title(['E_{mag} for m =' num2str(m)  ' vs \phi and \theta']);
    saveas(fig2,['./plots/space_time_2D_far_field/' figname '.png']); % Save the figure as a PNG file
end


%% Different ST Codes
function st_code = generate_special_st_code()
    % This function generates the space time code for the Finer Beam Control example in our report 
    % Initialize 8x8x8 matrix
    % Dimensions: 8(X) x 8(Y) x 8(Time)
    st_code = zeros(8, 8, 8);
    
    % Rows 1-8 correspond to Time-coding sequences 1-8
    pattern_2d = [
        2, 3, 0, 0, 1, 1, 2, 2; % T=1
        3, 2, 3, 0, 0, 1, 1, 2; % T=2
        3, 3, 0, 1, 1, 1, 2, 2; % T=3
        2, 0, 0, 0, 1, 1, 1, 2; % T=4
        2, 3, 0, 3, 0, 1, 2, 1; % T=5
        2, 3, 3, 0, 1, 2, 2, 2; % T=6
        3, 3, 3, 0, 0, 1, 1, 3; % T=7
        3, 3, 3, 0, 0, 0, 1, 2  % T=8
    ];
    
    % We replicate the 2D pattern for every X-element
    for x = 1:8
        st_code(x, :, :) = transpose(pattern_2d);
    end
end


