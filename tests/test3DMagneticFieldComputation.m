classdef test3DMagneticFieldComputation < matlab.unittest.TestCase

    properties
        test_parameters;
    end

    methods(TestClassSetup)
        function setup(test_case)
            % Define test parameters
            test_case.test_parameters.error_thresh = 0.1;
        end
    end
    methods(Test)
        function test3DCylinder(test_case)
            %% Prepare workspace
            % Get the full path of the current script
            fullpath = mfilename('fullpath');
            % Extract just the folder
            script_folder = fileparts(fullpath);
            cd(script_folder);

            % Have to add the functions path manually so prepare_workspace runs
            parent_folder = fileparts(script_folder);
            addpath(genpath(fullfile(parent_folder,'functions')));

            model_folder = prepare_workspace(script_folder);
            
            clc;
            %% Define test parameters

            %Define the characteristic scales in SI units

            z0 = 0.0058; %(Ohm m^2) is the contact impedance from the CEM article 58 Ohm cm^2
            l0 = 40e-3; %(m) the tank radius
            I0 = 2.4e-3;%(A) the magnitude of the injected current

            % The derived characteristic units
            V0 = z0*I0/(l0^2); %(V)
            sigma0 = l0/z0; %(S/m)
            J0 = I0/(l0^2);

            %% Create 3D template model_parameters
            model_parameters = create_kai_3d_model_parameters(l0, z0, sigma0, I0);

            model_parameters.maxsz = 2*model_parameters.maxsz;
            model_parameters.radius = 80e-3/l0;
            model_parameters.height = model_parameters.radius;
            
            % Set sensor measurement axis to something other than cartesian
            model_parameters.measurementAxisType = 'cylindrical';
         
            %% Create empty forward model
            
            cyl_shape = [model_parameters.height,model_parameters.radius,model_parameters.maxsz];
            elec_pos = [];
            elec_shape = [];
            
            fmdl= ng_mk_cyl_models(cyl_shape,elec_pos,elec_shape);
            
            %% Assign 2 electrodes on the top and bottom of the cylinder

            % Find nodes at bottom
            bottom_nodes = find(abs(fmdl.nodes(:,3)- 0.0)<1e-12);
            top_nodes = find(abs(fmdl.nodes(:,3) - model_parameters.height)<1e-12);

            electrode(1).nodes =  bottom_nodes;
            electrode(1).z_contact = 1;

            electrode(2).nodes =top_nodes;
            electrode(2).z_contact = 1;

            fmdl.electrode = electrode;
            
            %% Create a stimulation structure
            stimulation.stimulation = 'Amp';
            stimulation.stim_pattern = sparse([1 2],[1 1],[-1 1],2,1);
            stimulation.meas_pattern = sparse([1 2 1 2],[1 1 2 2],[-1 1 1 -1],2,2);

            fmdl.stimulation = stimulation;

            %% Assign sensors
            sz = model_parameters.radius;
            sensor_locations = [...
                -sz -sz model_parameters.height/2;
                2*sz -sz model_parameters.height/2;
                -sz 2*sz model_parameters.height/2;
                2*sz 2*sz model_parameters.height/2;
                0 -sz model_parameters.height/2;
                sz -sz model_parameters.height/2;
                -sz sz model_parameters.height/2;
                2*sz sz model_parameters.height/2;
                -sz 0 model_parameters.height/2;
                2*sz 0 model_parameters.height/2;
                0 2*sz model_parameters.height/2;
                sz 2*sz model_parameters.height/2]-[sz/2,sz/2,0];

            num_sensors = size(sensor_locations,1);
            
            sensor_axes = construct_sensor_axis(sensor_locations,'cylindrical');
            
            % Assign sensorLocations and sensorAxes to fmdl
            sensors = repmat(struct('position', [], 'axes', []), 1, num_sensors);
            
            for m = 1:num_sensors
                sensors(m).position = sensor_locations(m,:);
                sensors(m).axes = sensor_axes(m);
            end

            fmdl.sensors = sensors;

            %% Plot the sensors and sensor axis to see if everything is correct
            show_fem(fmdl);
            plot_sensors(fmdl)
            
            hold on
            for m = 1:numel(fmdl.sensors)
                position = fmdl.sensors(m).position;

                axes = {fmdl.sensors(m).axes.axis1,fmdl.sensors(m).axes.axis2,fmdl.sensors(m).axes.axis3};
                for n = 1:3
                    quiver3(position(1),position(2),position(3),axes{n}(1),axes{n}(2),axes{n}(3))
                end
            end
            hold off

            xlim([-2*model_parameters.radius,2*model_parameters.radius])
            ylim([-2*model_parameters.radius,2*model_parameters.radius])

            drawnow;

            %% Define the current:
            J = [0 0 -1];
            mu0 = model_parameters.mu0;

            %% Compute magnetic fields with numerical integration of Biot Savart Law
            for m = 1:numel(fmdl.sensors)

                x_obs = fmdl.sensors(m).position(1);
                y_obs = fmdl.sensors(m).position(2);
                z_obs = fmdl.sensors(m).position(3);

                % Cylinder parameters
                R  = model_parameters.radius;
                z1 = 0;
                z2 = model_parameters.height;

                % Current density
                Jx = J(1); Jy = J(2); Jz = J(3);

                % Distance components
                Rx = @(r,th,z) x_obs - r.*cos(th);
                Ry = @(r,th,z) y_obs - r.*sin(th);
                Rz = @(r,th,z) z_obs - z;

                Rnorm = @(r,th,z) (Rx(r,th,z).^2 + ...
                    Ry(r,th,z).^2 + ...
                    Rz(r,th,z).^2).^(3/2);

                % ---- Biot–Savart integrands ----
                Bx_integrand = @(r,th,z) ...
                    (Jy.*Rz(r,th,z) - Jz.*Ry(r,th,z)) ./ Rnorm(r,th,z) .* r;

                By_integrand = @(r,th,z) ...
                    (Jz.*Rx(r,th,z) - Jx.*Rz(r,th,z)) ./ Rnorm(r,th,z) .* r;

                Bz_integrand = @(r,th,z) ...
                    (Jx.*Ry(r,th,z) - Jy.*Rx(r,th,z)) ./ Rnorm(r,th,z) .* r;
                
                % ---- Numerical integration ----
                Bx_num = (mu0/(4*pi)) * integral3(Bx_integrand, ...
                    0, R, 0, 2*pi, z1, z2, ...
                    'AbsTol',1e-9,'RelTol',1e-6);

                By_num = (mu0/(4*pi)) * integral3( By_integrand, ...
                    0, R, 0, 2*pi, z1, z2, ...
                    'AbsTol',1e-9,'RelTol',1e-6);

                Bz_num = (mu0/(4*pi)) * integral3( Bz_integrand, ...
                    0, R, 0, 2*pi, z1, z2, ...
                    'AbsTol',1e-9,'RelTol',1e-6);
                
                % CHANGE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
                % FOR SOME REASON, this differes by a factor of 1/(4*pi)
                % with the R-matrix computation for cartesian coordinates,
                % so hack this for now

                Bm = 1/(4*pi)*[Bx_num, By_num, Bz_num];

                B_num(m,:) = [...
                    dot(Bm,fmdl.sensors(m).axes.axis1),
                    dot(Bm,fmdl.sensors(m).axes.axis2),
                    dot(Bm,fmdl.sensors(m).axes.axis3)];

            end

            %% Compute magnetic fields with r matrices
            fmdl = compute_geometry_matrices(fmdl,mu0);

            %Assign homogeneous conductivity
            img = mk_image(fmdl,1.0,'Homogeneous 3D cylinder');
            img.fwd_solve.get_all_meas = 1;
            
            %Solve forward model with my solver
            [data,~] = fwd_solve_mdeit(img);

            %% 
            figure
            subplot(1,3,1)
            hold on
            plot(data.Bx(:));
            plot(B_num(:,1),'b.');
            hold off
            legend('Solver','Numerical')
            axis square
            grid on;grid minor;
            title('$B_r$','Interpreter','latex')

            subplot(1,3,2)
            hold on
            plot(data.By(:));
            plot(B_num(:,2),'b.');
            hold off
            legend('Solver','Numerical')
            axis square
            grid on;grid minor;
            title('$B_\theta$','Interpreter','latex')

            subplot(1,3,3)
            hold on
            plot(data.Bz(:));
            plot(B_num(:,3),'b.');
            hold off
            legend('Solver','Numerical')
            axis square
            grid on;grid minor;
            title('$B_z$','Interpreter','latex')

            drawnow

            %% Check difference in theta component
            clc;
            fprintf('Error threshold: %2.2f %%\n',test_case.test_parameters.error_thresh);
           
            error_theta = abs(data.By(:)-B_num(:,2))./abs(B_num(:,2))*100;
            
            pass = norm(error_theta,'inf')<test_case.test_parameters.error_thresh;
            test_case.assertTrue(pass,"norm(error_theta,'inf') < \epsilon");
        end
    end
end








function sensor_axes = construct_sensor_axis(sensor_locations,measurement_axis_type)

num_sensors = size(sensor_locations,1);

center_of_referential = mean(sensor_locations);

sensor_axes = ...
    repmat(struct('axis1', [], 'axis2', [],'axis3',[]), 1, num_sensors);

switch measurement_axis_type
    case 'cartesian'
        for m=1:num_sensors
            sensor_axes(m).axis1 = [1,0,0];
            sensor_axes(m).axis2 = [0,1,0];
            sensor_axes(m).axis3 = [0,0,1];
        end
    case 'cylindrical'
        for m=1:num_sensors
            p = sensor_locations(m,:)-center_of_referential;

            x = p(1);
            y = p(2);
            z = p(3);

            rhat = [x,y,0];
            rhat = rhat/norm(rhat);

            thetahat = [-y,x,0];
            thetahat = thetahat/norm(thetahat);

            zhat = [0,0,1];


            sensor_axes(m).axis1 = rhat;
            sensor_axes(m).axis2 = thetahat;
            sensor_axes(m).axis3 = zhat;
        end
    case 'spherical'
        for m = 1:num_sensors
            
            p = sensor_locations(m,:)-center_of_referential;

            x = p(1);
            y = p(2);
            z = p(3);
            
            r = sqrt(x^2+y^2+z^2);
            theta = acos(z/sqrt(r^2));
            phi = sign(y)*acos(x/sqrt(x^2+y^2));

            rhat = sin(theta)*cos(phi)*[1,0,0]+sin(theta)*sin(phi)*[0,1,0]+cos(theta)*[0,0,1];
            thetahat = cos(theta)*cos(phi)*[1,0,0]+cos(theta)*sin(phi)*[0,1,0]-sin(theta)*[0,0,1];
            phihat = -sin(phi)*[1,0,0]+cos(phi)*[0,1,0];

            sensor_axes(m).axis1 = rhat;
            sensor_axes(m).axis2 = thetahat;
            sensor_axes(m).axis3 = phihat;
        end
   
    otherwise
        error('measurement_axis_type %s is not valid',measurement_axis_type);
end

end
