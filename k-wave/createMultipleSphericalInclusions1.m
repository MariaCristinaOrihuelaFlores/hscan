function matrix = createMultipleSphericalInclusions1(Nx, Ny, Nz, centers, radii, numSpheres)
    % N is the size of the NxNxN matrix
    % centers is an Mx3 matrix, where each row represents the x, y, z center of a sphere
    % radii is a vector of length M, containing the radius of each sphere
    % numSpheres is the number of spheres
    
    % Initialize the 3D matrix with zeros
    matrix = zeros(Nx, Ny, Nz);
    
    % Loop through each sphere
    for s = 1:numSpheres
        % Extract the center for the current sphere
        x0 = centers(s, 1);
        y0 = centers(s, 2);
        z0 = centers(s, 3);
        
        % Calculate the distance from each point in the grid to the center of the current sphere
        for i = 1:Nx
            for j = 1:Ny
                for k = 1:Nz
                    % Calculate the distance from the current point to the center of the sphere
                    dist = sqrt((i - x0)^2 + (j - y0)^2 + (k - z0)^2);
                    
                    % Check if the distance is within the radius of the sphere
                    if dist <= radii(s)
                        % If it is, set the corresponding point in the matrix to 1
                        matrix(i, j, k) = 1;
                    end
                end
            end
        end
    end
end
