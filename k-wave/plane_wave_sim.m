%% Plane wave


clear; clc; close all;
% simulation settings
DATA_CAST       = 'single';     % set to 'single' or 'gpuArray-single' to speed up computations
RUN_SIMULATION  = true;         % set to false to reload previous results instead of running simulation 
% =========================================================================
% DEFINE THE K-WAVE GRID
% =========================================================================

% set the size of the perfectly matched layer (PML)
pml_x_size = 40;                % [grid points]
pml_y_size = 10;                % [grid points]
pml_z_size = 10;                % [grid points]

% set total number of grid pointsta not including the PML
Nx = 1248 - 2 * pml_x_size;      % [grid points]
Ny = 675 - 2 * pml_y_size;      % [grid points]
Nz = 675 - 2 * pml_z_size;      % [grid points]

% set desired grid size in the x-direction not including the PML
x = 70.08e-3;                      % [m]

% calculate the spacing between the grid points
dx = x / Nx;                    % [m]
dy = dx;                        % [m]
dz = dx;                        % [m]

% create the k-space grid
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% =========================================================================
% DEFINE THE MEDIUM PARAMETERS
% =========================================================================

% define the properties of the propagation medium
c0 = 1540;                      % [m/s]
rho0 = 1000;                    % [kg/m^3]
medium.alpha_coeff = 0;      % [dB/(MHz^y cm)]
medium.alpha_power = 0;

% create the time array
t_end = (Nx * dx) * 2.2 / c0;   % [s]
kgrid.makeTime(c0, [], t_end);
 
% =========================================================================
% DEFINE THE INPUT SIGNAL
% =========================================================================

% define properties of the input signal
source_strength = 1e6;          % [Pa]
tone_burst_freq = 7.6e6;        % [Hz]
tone_burst_cycles = 4;

% create the input signal using toneBurst 
input_signal = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);

% scale the source magnitude by the source_strength divided by the
% impedance (the source is assigned to the particle velocity)
input_signal = (source_strength ./ (c0 * rho0)) .* input_signal;
 
% =========================================================================
% DEFINE THE ULTRASOUND TRANSDUCER
% =========================================================================

% physical properties of the transducer
transducer.number_elements = 128;  	% total number of transducer elements
transducer.element_width = round(270e-6/dx);       % width of each element [grid points]
transducer.element_length = Nz-2;  	% length of each element [grid points]
transducer.element_spacing = 0;  	% spacing (kerf  width) between the elements [grid points]
transducer.radius = inf;            % radius of curvature of the transducer [m]

% calculate the width of the transducer in grid points
transducer_width = transducer.number_elements * transducer.element_width ...
    + (transducer.number_elements - 1) * transducer.element_spacing;

% use this to position the transducer in the middle of the computational grid
transducer.position = round([1, Ny/2 - transducer_width/2, Nz/2 - transducer.element_length/2]);

% properties used to derive the beamforming delays
transducer.sound_speed = c0;                    % sound speed [m/s]
transducer.focus_distance = inf;              % focus distance [m]
transducer.elevation_focus_distance = inf;    % focus distance in the elevation plane [m]
transducer.steering_angle = 0;                  % steering angle [degrees]

% apodization
transducer.transmit_apodization = 'Hanning';    
transducer.receive_apodization = 'Rectangular';

% define the transducer elements that are currently active
number_active_elements = 8;
transducer.active_elements = ones(transducer.number_elements, 1);

% append input signal used to drive the transducer
transducer.input_signal = input_signal;

% create the transducer using the defined settings
transducer = kWaveTransducer(kgrid, transducer);

% print out transducer properties
transducer.properties; 

% =========================================================================
% DEFINE THE MEDIUM PROPERTIES
% =========================================================================

Nx_tot = Nx;
Ny_tot = Ny;
Nz_tot = Nz;

% define a random distribution of scatterers for the medium
background_map_mean = 1;
background_map_std = 0.008;
background_map = background_map_mean + background_map_std * randn([Nx_tot, Ny_tot, Nz_tot]);

% define a random distribution of scatterers for the highly scattering
%region
scattering_map = randn([Nx_tot, Ny_tot, Nz_tot]);
scattering_c0 = c0 + 410 + 75 * scattering_map;
scattering_c0(scattering_c0 > 2010) = 2010;
scattering_c0(scattering_c0 < 1810) = 1810;
scattering_rho0 = scattering_c0 / 2.03;

% define properties
sound_speed_map = c0 * ones(Nx_tot, Ny_tot, Nz_tot) .* background_map;
density_map = rho0 * ones(Nx_tot, Ny_tot, Nz_tot) .* background_map;

% define a sphere for a highly scattering region
radius = (5e-3)/2;      % [m]
x_pos = 23.36e-3;    % [m]
y_pos = 19.65e-3;    % [m]
scattering_region1 = createMultipleSphericalInclusions1(Nx_tot, Ny_tot, Nz_tot, [round(x_pos/dx);round(y_pos/dx);Nz_tot/2]', round(radius/dx), 1);
%scattering_region1 = makeBall(Nx_tot, Ny_tot, Nz_tot, round(x_pos/dx), round(y_pos/dx), Nz_tot/2, round(radius/dx)); 

% assign region
sound_speed_map(scattering_region1 == 1) = scattering_c0(scattering_region1 == 1);
density_map(scattering_region1 == 1) = scattering_rho0(scattering_region1 == 1);

% define a sphere for a highly scattering region
radius = (10e-3)/2;      % [m]
x_pos = 46.72e-3;    % [m]
y_pos = 19.65e-3;      % [m]
scattering_region2 = createMultipleSphericalInclusions1(Nx_tot, Ny_tot, Nz_tot, [round(x_pos/dx);round(y_pos/dx);Nz_tot/2]', round(radius/dx), 1);

% assign region
sound_speed_map(scattering_region2 == 1) = scattering_c0(scattering_region2 == 1);
density_map(scattering_region2 == 1) = scattering_rho0(scattering_region2 == 1); 
medium.sound_speed = sound_speed_map(:,:, :);
medium.density = density_map(:, :, :);
%% 

% =========================================================================
% RUN THE SIMULATION
% =========================================================================
input_args = {...
    'PMLInside', false, 'PMLSize', [pml_x_size, pml_y_size, pml_z_size], ...
    'DataCast', DATA_CAST, 'DataRecast', true, 'PlotSim', false};

% run the simulation
[sensor_data] = kspaceFirstOrder3D(kgrid, medium, transducer, transducer, input_args{:});
scan_lines = transducer.combine_sensor_data(sensor_data) ;
save us_bmode_plane_waves sensor_data 
    
%% 
% =========================================================================
% PROCESS THE RESULTS
% =========================================================================

% -----------------------------
% Remove Input Signal
% -----------------------------

% create a window to set the first part of each scan line to zero to remove
% interference from the input signal
scan_line_win = getWin(kgrid.Nt * 2, 'Tukey', 'Param', 0.05).';
scan_line_win = [zeros(1, length(input_signal) * 2), scan_line_win(1:end/2 - length(input_signal) * 2)];

% apply the window to each of the scan lines
scan_lines = bsxfun(@times, scan_line_win, scan_lines);

% store a copy of the middle scan line to illustrate the effects of each
% processing step
scan_line_example(1, :) = scan_lines(end/2, :);

% -----------------------------
% Time Gain Compensation
% -----------------------------

% create radius variable assuming that t0 corresponds to the middle of the
% input signal
t0 = length(input_signal) * kgrid.dt / 2;
r = c0 * ( (1:length(kgrid.t_array)) * kgrid.dt - t0 ) / 2;    % [m]

% define absorption value and convert to correct units
tgc_alpha_db_cm = medium.alpha_coeff * (tone_burst_freq * 1e-6)^medium.alpha_power;
tgc_alpha_np_m = tgc_alpha_db_cm / 8.686 * 100;

% create time gain compensation function based on attenuation value and
% round trip distance
tgc = exp(tgc_alpha_np_m * 2 * r);

% apply the time gain compensation to each of the scan lines
scan_lines = bsxfun(@times, tgc, scan_lines);

% store a copy of the middle scan line to illustrate the effects of each
% processing step
scan_line_example(2, :) = scan_lines(end/2, :);

% -----------------------------
% Frequency Filtering
% -----------------------------

% filter the scan lines using both the transmit frequency and the second
% harmonic
scan_lines_fund = gaussianFilter(scan_lines, 1/kgrid.dt, tone_burst_freq, 100, true);
set(gca, 'XLim', [0, 6]);
scan_lines_harm = gaussianFilter(scan_lines, 1/kgrid.dt, 2 * tone_burst_freq, 30, true);
set(gca, 'XLim', [0, 6]);

% store a copy of the middle scan line to illustrate the effects of each
% processing step
scan_line_example(3, :) = scan_lines_fund(end/2, :);

% -----------------------------
% Envelope Detection
% -----------------------------

% envelope detection
scan_lines_fund = envelopeDetection(scan_lines_fund);
scan_lines_harm = envelopeDetection(scan_lines_harm);

% store a copy of the middle scan line to illustrate the effects of each
% processing step
scan_line_example(4, :) = scan_lines_fund(end/2, :);

% -----------------------------
% Log Compression
% -----------------------------

% normalised log compression
compression_ratio = 3;
scan_lines_fund = logCompression(scan_lines_fund, compression_ratio, true);
scan_lines_harm = logCompression(scan_lines_harm, compression_ratio, true);

% store a copy of the middle scan line to illustrate the effects of each
% processing step
scan_line_example(5, :) = scan_lines_fund(end/2, :);

% -----------------------------
% Scan Conversion
% -----------------------------

% upsample the image using linear interpolation
scale_factor = 2;
scan_lines_fund = interp2(1:kgrid.Nt, (1:number_scan_lines).', scan_lines_fund, 1:kgrid.Nt, (1:1/scale_factor:number_scan_lines).');
scan_lines_harm = interp2(1:kgrid.Nt, (1:number_scan_lines).', scan_lines_harm, 1:kgrid.Nt, (1:1/scale_factor:number_scan_lines).');
%% 
% =========================================================================
% VISUALISATION OF SIMULATION LAYOUT
% =========================================================================

% % uncomment to generate a voxel plot of the simulation layout
% 
% % physical properties of the transducer
% transducer_plot.number_elements = 32 + number_scan_lines - 1;
% transducer_plot.element_width = 2;
% transducer_plot.element_length = 24;
% transducer_plot.element_spacing = 0;
% transducer_plot.radius = inf;
% 
% % transducer position
% transducer_plot.position = round([1, Ny/2 - transducer_width/2, Nz/2 - transducer.element_length/2]);
% 
% % create expanded grid
% kgrid_plot = kWaveGrid(Nx_tot, dx, Ny_tot, dy, Nz, dz);
% kgrid_plot.setTime(1, 1);
% 
% % create the transducer using the defined settings
% transducer_plot = kWaveTransducer(kgrid_plot, transducer_plot);
% 
% % create voxel plot of transducer mask and 
% voxelPlot(single(transducer_plot.active_elements_mask | scattering_region1 | scattering_region2 | scattering_region3));
% view(26, 48);

