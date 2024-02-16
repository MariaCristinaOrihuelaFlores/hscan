 define properties
sound_speed_map = c0_min * ones(Nx_tot, Ny_tot, Nz_tot) .* background_map;
density_map = rho0 * ones(Nx_tot, Ny_tot, Nz_tot) .* background_map;
medium.sound_speed = sound_speed_map(:,:, :);
medium.density = density_map(:, :, :);
% % =========================================================================
% % RUN THE SIMULATION
% % =========================================================================

% set the input settings
input_args = {...
'PMLInside', false, 'PMLSize', [pml_x_size, pml_y_size, pml_z_size], ...
'PMLAlpha', [pml_x_alpha, pml_y_alpha, pml_z_alpha], 'DataCast', DATA_CAST, 'DataRecast', true, 'PlotSim', false};
% save input files to disk

% run the simulation
[sensor_data] = kspaceFirstOrder3DG(kgrid, medium, transducer, transducer, input_args{:});
scan_lines = transducer.combine_sensor_data(sensor_data) ;
%
save sensor_data
%%
figure;imagesc(sound_speed_map(:,:,round(Nz/2)))