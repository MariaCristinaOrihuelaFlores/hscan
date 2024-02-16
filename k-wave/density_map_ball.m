function sum_ball = density_map_ball (Nx,Ny,Nz)
for m=1:1 %concentrations
   %concentration = 0.002 + (0.01 - 0.002) * rand(1, 1); % [%]
   for c = 1:1
       diameter = 20 + (8 - 6) * rand(1, 1); %[grid]
       radii=diameter/2;
       num_points = 5; %round(concentration * Nx * Ny * Nz)
       cx = randi(Nx, 1, num_points);
       cy = randi(Ny, 1, num_points);
       cz = randi(Nz, 1, num_points);
       sum_ball = createMultipleSphericalInclusions1(Nx, Ny, Nz, [cx;cy;cz]', radii.*ones(1,num_points), num_points);
       %figure;voxelPlot(sum_ball);
   end
end
