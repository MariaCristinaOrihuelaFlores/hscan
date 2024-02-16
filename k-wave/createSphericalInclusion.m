function matrix = createSphericalInclusion(Nx,Ny,Nz, cx, cy, cz, r)
    % Initialize the 3D matrix with zeros
    matrix = zeros(Nx, Ny, Nz);
    num_points = length(cx);
    for num = 1:num_points
        % Loop through each element in the matrix
        for x = 1:Nx
            for y = 1:Ny
                for z = 1:Nz
                    % Calculate the distance from the center
                    dist = sqrt((x - x0)^2 + (y - y0)^2 + (z - z0)^2);
                    
                    % Check if the distance is less than or equal to the radius
                    if dist <= r
                        % If yes, then it's inside the sphere, so set the element to 1
                        matrix(x, y, z) = 1;
                    end
                end
            end
        end
    end
end