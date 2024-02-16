function matrix = createMultipleSphericalInclusions(Nx,Ny,Nz, centers, radii)
    % N is the size of the NxNxN matrix
    % centers is an Mx3 matrix, where each row represents the x, y, z center of a sphere
    % radii is a vector of length M, containing the radius of each sphere
    
    % Initialize the 3D matrix with zeros
    matrix = zeros(Nx, Ny, Nz);
    
    % Generate a 3D grid of coordinates
    [X, Y, Z] = meshgrid(1:Nx, 1:Ny, 1:Nz);
    
    % Loop through each sphere
    for i = 1:length(radii)
        % Extract the center for the current sphere
        x0 = centers(i, 1);
        y0 = centers(i, 2);
        z0 = centers(i, 3);
        
        % Calculate the distance from each point in the grid to the center of the current sphere
        dist = sqrt((X - x0).^2 + (Y - y0).^2 + (Z - z0).^2);
        
        % Find the points where the distance is less than or equal to the radius and set them to 1
        matrix(dist <= radii(i)) = 1;
    end
end