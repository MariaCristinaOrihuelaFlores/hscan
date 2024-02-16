clear; clc; close all;
% simulation settings
DATA_CAST       = 'single';     % set to 'single' or 'gpuArray-single' to speed up computations
RUN_SIMULATION  = true;         % set to false to reload previous results instead of running simulation
% Grid
% =========================================================================
% DEFINE THE K-WAVE GRID
% =========================================================================
    
% set the size of the perfectly matched layer (PML)
pml_x_size = 40;                % [grid points] *2
pml_y_size = 10;                % [grid points]
pml_z_size = 10;                % [grid points]
            
Nx = 400 - 2 * pml_x_size;      % [grid points] 2048
Ny = 400- 2 * pml_y_size;      % [grid points] 1024
Nz = 100 - 2 * pml_z_size;      % [grid points]
            
% calculate the spacing between the grid points
ratio = 1; % grid per element [grid points] ?
pitch =  0.3e-3;
dx = pitch/ratio;           % [m]
dy = dx;
dz = dx;
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);
%
% =========================================================================
% DEFINE THE MEDIUM PARAMETERS
% =========================================================================
alpha = 0;
cycles= 4.5;
% define the properties of the propagation medium
c0 = 1540;                      % [m/s]
rho0 = 1000;                    % [kg/m^3]
medium.alpha_coeff = alpha;      % [dB/(MHz^y cm)]
medium.alpha_power = 0;


% create the time array
t_end = (Nx * dx) * 2.2 / c0;   % [s]
kgrid.makeTime(c0, [], t_end);


%
% =========================================================================
% DEFINE THE INPUT SIGNAL
% ========================================================================
% define properties of the input signal
source_strength = 1e6;          % [Pa]
tone_burst_freq = 5e6;        % [Hz]
tone_burst_cycles = cycles; %5 10 15 20  

input_signal_norm = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles, 'Envelope', 'Gaussian');
input_signal = (source_strength ./ (c0 * rho0)) .* input_signal_norm;
%% 
% =========================================================================
% DEFINE THE ULTRASOUND TRANSDUCER
% =========================================================================
        
% physical properties of the transducer
transducer.number_elements = 128;  	% total number of transducer elements
transducer.element_width = floor(Ny/transducer.number_elements);       % width of each element [grid points] 9
transducer.element_length = Nz-1;  	% length of each element [grid points] 166
transducer.element_spacing = 0;  	% spacing (kerf  width) between the elements [grid points] 1
transducer.radius = inf;            % radius of curvature of the transducer [m]
        
% calculate the width of the transducer in grid points
transducer_width = transducer.number_elements * transducer.element_width ...
            + (transducer.number_elements - 1) * transducer.element_spacing;
        
% use this to position the transducer in the middle of the computational grid
transducer.position = round([1, Ny/2 - transducer_width/2, Nz/2 - transducer.element_length/2]);
        
% properties used to derive the beamforming delays
transducer.sound_speed = c0;                    % sound speed [m/s]
transducer.focus_distance = inf;              % focus distance [m]
transducer.elevation_focus_distance = inf;    % focus distance in the elevation plane [m] 0.02
transducer.steering_angle = 0;                  % steering angle [degrees]
        
% apodization
%transducer.transmit_apodization = 'Rectangular';    
%transducer.receive_apodization = 'Rectangular';
        
% define the transducer elements that are currently active
number_active_elements = 8;
transducer.active_elements = ones(transducer.number_elements, 1);
        
% append input signal used to drive the transducer
transducer.input_signal = input_signal;
        
% create the transducer using the defined settings
transducer = kWaveTransducer(kgrid, transducer);

%transducer.plot
% print out transducer properties
%transducer.properties;
%% 
% =========================================================================
% DEFINE THE MEDIUM PROPERTIES
% =========================================================================
        
% define a large image size to move across
% number_scan_lines = 500;
Nx_tot = Nx;
% Ny_tot = Ny + number_scan_lines * transducer.element_width;
Ny_tot=Ny;
Nz_tot = Nz;
        
alpha_coeff_map = alpha; 

mask = density_map_ball(Nx,Ny,Nz);
sos = zeros(Nx, Ny, Nz);
dens = zeros(Nx, Ny, Nz);
sos(mask==1) = 1950; sos(mask==0) = 1540;
dens(mask==1) = 960;dens(mask==0) = 1000;

medium.sound_speed = sos;
%% 
% =========================================================================
% RUN THE SIMULATION
% =========================================================================
input_args = {...
'PMLInside', false, 'PMLSize', [pml_x_size, pml_y_size, pml_z_size], ...
'DataCast', DATA_CAST, 'DataRecast', true, 'PlotSim', false};
        
% run the simulation if set to true, otherwise, load previous results from
% disk
% if RUN_SIMULATION
        
% set medium position
%medium_position = 1;
number_scan_lines = Ny;
medium.density = dens;
scan_lines = [];

medium.density = dens;
sensor_data = kspaceFirstOrder3D(kgrid, medium, transducer, transducer, input_args{:});
% for scan_line_index = 1:number_scan_lines% number slices
%     disp('');
%     disp(['Computing scan line PLANE WAVE...']);
%     medium.density = dens(:, medium_position:medium_position + Ny - 1, :);
%     % run the simulation
%     sensor_data = kspaceFirstOrder3D(kgrid, medium, transducer, transducer, input_args{:});
%     
%     % extract the scan line from the sensor data
%     scan_lines(scan_line_index, :) = transducer.scan_line(sensor_data);
%         
%     % update medium position
%     medium_position = medium_position + transducer.element_width;
% end
%% das
RF = sensor_data;
%% 

file_out = 'RF_prueba2';
save(file_out,'RF'); %clear transducer
%% 
scan_lines_fund = envelopeDetection(RF);
%% 
scale_factor = 2; 
scan_lines_fund1 = interp2(1:kgrid.Nt, (1:128).', scan_lines_fund, 1:kgrid.Nt, (1:1/scale_factor:number_scan_lines).'); 

%% 
figure;imagesc(scan_lines_fund1);colormap gray;