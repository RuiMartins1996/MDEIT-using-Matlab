classdef test2DMagneticFieldComputation < matlab.unittest.TestCase


    properties
        test_parameters;
    end

    methods(TestClassSetup)
        function setup(test_case)
            % Define test parameters
            test_case.test_parameters.mu0 = 4*pi*1e-7;
        end
    end

    methods(TestMethodSetup)
        % Setup for each test
    end

    methods(Test)
        % Test methods
        function test2DConductivePlate(test_case)
            %% Prepare workspace
            % Get the full path of the current script
            fullpath = mfilename('fullpath');
            % Extract just the folder
            script_folder = fileparts(fullpath);
            cd(script_folder);

            % Have to add the functions path manually so prepare_workspace runs
            parent_folder =fileparts(script_folder);
            addpath(genpath(fullfile(parent_folder,'functions')));

            prepare_workspace(script_folder);

            %% Create a square conductive plaque
            sz = 1.0;           % size
            loop = [0 0;  sz 0; sz sz; 0 sz];
            maxsz = 0.05;       % element size

            fmdl = create_2d_conductive_square(loop,maxsz);

            show_fem(fmdl);
            %% Add a single stimulation pattern
            current_amplitude = 1;

            num_electrodes = 2;
            num_meas = 1;

            stimulation.stimulation = 'Amp';
            % Current flowing from bottom to top: J = (0,+Jy,0);
            stimulation.stim_pattern = sparse(num_electrodes,1);
            stimulation.stim_pattern(1,1) = -current_amplitude;
            stimulation.stim_pattern(2,1) = +current_amplitude;

            stimulation.meas_pattern = sparse(num_meas,num_electrodes);
            stimulation.meas_pattern(1,1) = 1;
            stimulation.meas_pattern(1,1) = -1;

            fmdl.stimulation = stimulation;
            %% Complete model with mdeit structures
            mu0 = test_case.test_parameters.mu0;

            sensor_positions = [...
                -sz -sz 0;2*sz -sz 0;-sz 2*sz 0; 2*sz 2*sz 0;...
                0 -sz 0;sz -sz 0;-sz sz 0; 2*sz sz 0;...
                -sz 0 0;2*sz 0 0;0 2*sz 0; sz 2*sz 0];

            fmdl = assign_magnetometers(fmdl,sensor_positions);
            fmdl = compute_geometry_matrices(fmdl,mu0);

            plot_sensors(fmdl)
            axis([-2*sz,3*sz,-2*sz,3*sz])

            %% Assign homogeneous conductivity
            img = mk_image(fmdl,1.0,'Homogeneous 2D square');
            img.fwd_solve.get_all_meas = 1;
            %% Solve forward model with my solver
            [data,u] = fwd_solve_mdeit(img);
            %%
            Jy = -1.0; %for this solution d_y u = -1, so Jy = -\sigma d_y u = 1.0; SO WHY -1?
            J = [0,Jy,0];
            %% Solve forward problem with numerical integration of Biot Savart law
            for m = 1:size(sensor_positions,1)

                x_obs = sensor_positions(m,1);
                y_obs = sensor_positions(m,2);
                z_obs = sensor_positions(m,3);

                % ---- Numerical integration of Biot–Savart ----
                x1 = loop(1); x2 = loop(2);
                y1 = loop(3); y2 = loop(4);
                Jx = J(1); Jy = J(2);

                % Define integrands for each component of B
                Bx_integrand = @(xp, yp) Jy .* (z_obs) ./ ((x_obs - xp).^2 + (y_obs - yp).^2 + z_obs^2).^(3/2);
                By_integrand = @(xp, yp) -Jx .* (z_obs) ./ ((x_obs - xp).^2 + (y_obs - yp).^2 + z_obs^2).^(3/2);
                Bz_integrand = @(xp, yp) (Jx .* (y_obs - yp) - Jy .* (x_obs - xp)) ./ ...
                    ((x_obs - xp).^2 + (y_obs - yp).^2 + z_obs^2).^(3/2);

                % Perform numerical integrations
                Bx_num = (mu0/(4*pi)) * integral2(Bx_integrand, x1, x2, y1, y2, ...
                    'AbsTol',1e-9, 'RelTol',1e-6);
                By_num = (mu0/(4*pi)) * integral2(By_integrand, x1, x2, y1, y2, ...
                    'AbsTol',1e-9, 'RelTol',1e-6);
                Bz_num = (mu0/(4*pi)) * integral2(Bz_integrand, x1, x2, y1, y2, ...
                    'AbsTol',1e-9, 'RelTol',1e-6);

                B_num(m,:) = [Bx_num, By_num, Bz_num];
            end

            %% Check if results are the same
            % normX = norm(data.Bx(:)-B_num(:,1),2);
            % normY = norm(data.By(:)-B_num(:,2),2);
            normZ = norm(data.Bz(:)-B_num(:,3),2);
            
            % Only the the z-direction component should be non-zero, so
            % check only that
            passZ = normZ/norm(B_num(:,3))*100 < 0.0001;

            test_case.assertTrue(passZ);
        end
    end
end


%% Functions
function fmdl = create_2d_conductive_square(loop,maxsz)

sz1 = max(loop(:,1))-min(loop(:,1));
sz2 = max(loop(:,2))-min(loop(:,2));

% Sanity check to see if its a square
if abs(sz1-sz2)>1e-3
    error('Expected a square');
else
    sz = sz1;
end

electrode_pos = [0.5 1;0.5,0];

electrode_witdh = 0.1;
electrode_maxsz = maxsz;
electrode_shape = [electrode_witdh,electrode_maxsz];

fmdl = ng_mk_2d_model({loop,maxsz},electrode_pos,electrode_shape);

% For some reason, NETGEN can't generate electrodes that
% intersect the square's corners, so lets add those manually

%Find the bottom and top nodes
bottom_nodes = find(abs(fmdl.nodes(:,2)-0)<maxsz/20);
top_nodes = find(abs(fmdl.nodes(:,2)-sz)<maxsz/20);

fmdl.electrode(1).nodes = top_nodes;
fmdl.electrode(2).nodes = bottom_nodes;
end


function fmdl = assign_magnetometers(fmdl,sensor_locations)

num_sensors = size(sensor_locations,1);

for m = 1:num_sensors
    sensor_axes(m).axis1 = [1,0,0];
    sensor_axes(m).axis2 = [0,1,0];
    sensor_axes(m).axis3 = [0,0,1];
end

% Assign sensorLocations and sensorAxes to fmdl
sensors = repmat(struct('position', [], 'axes', []), 1, num_sensors);

for m = 1:num_sensors
    sensors(m).position = sensor_locations(m,:);
    sensors(m).axes = sensor_axes(m);
end

fmdl.sensors = sensors;
end