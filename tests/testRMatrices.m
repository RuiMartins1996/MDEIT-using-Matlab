classdef testRMatrices < matlab.unittest.TestCase

    properties
        testParameters;
    end

    methods(TestClassSetup)
        % Shared setup for the entire test class
    end

    methods(TestMethodSetup)
        function setup(testCase)
            %% Define test parameters

            % Current density vector
            testCase.testParameters.J = [1 1 1];

            % Error threshold
            testCase.testParameters.errorThresh = 1;
        end
    end

    %% Test methods
    methods(Test)
        function testExactExpressionForRMatrices(testCase)
            %% Prepare workspace
            % Get the full path of the current script
            fullpath = mfilename('fullpath');
            % Extract just the folder
            script_folder = fileparts(fullpath);
            cd(script_folder);

            % Have to add the functions path manually so prepare_workspace runs
            parent_folder =fileparts(script_folder);
            addpath(genpath(fullfile(parent_folder,'functions')));

            model_folder = prepare_workspace(script_folder);

            rng(1);

            %% Extract Parameters
            errorThresh = testCase.testParameters.errorThresh;

            %% Model parameters
            z0 = 0.0058; %(Ohm m^2) is the contact impedance from the CEM article 58 Ohm cm^2
            l0 = 40e-3; %(m) the tank radius
            I0 = 2.4e-3;%(A) the magnitude of the injected current

            % The derived characteristic units
            V0 = z0*I0/(l0^2); %(V)
            sigma0 = l0/z0; %(S/m)
            J0 = I0/(l0^2);

            model_parameters = create_default_3d_model_parameters(l0, z0, sigma0, I0);

            model_parameters.maxsz = max(model_parameters.height,model_parameters.radius)*0.05;
            model_parameters.numOfElectrodesPerRing = 4;
            model_parameters.numOfSensors = 4;
            model_parameters.sensorRadius = model_parameters.radius*2;

            %% Make forward model

            [~,fmdls] = mk_mdeit_model(model_parameters,model_folder);

            fmdl = fmdls{1};

            sensorLocations = zeros(numel(fmdl.sensors),3);
            for i = 1: numel(fmdl.sensors)
                sensorLocations(i,:) = fmdl.sensors(i).position;
            end

            fprintf("Test Parameters: \n"+ ...
                "Error threshold: %g (%%)\n",testCase.testParameters.errorThresh);

            %% Compute R matrices with numerical method
            [RxNewtonCotes4,RyNewtonCotes4,RzNewtonCotes4] = testRMatrices.computeRmatricesNewtonCotes4(fmdl,sensorLocations);

            %% Compute R matrices with analytic method
            numElements = size(fmdl.elems,1);
            numSensors = size(sensorLocations,1);

            Rx = zeros(numSensors,numElements);
            Ry = zeros(numSensors,numElements);
            Rz = zeros(numSensors,numElements);

            figure;
            hold on;
            show_fem(fmdl);
            plot_sensors(fmdl);
            hold on

            countx = 0;
            county = 0;
            countz = 0;

            % It seems that elements 3737,3741 are the problem for these
            % situation

            elements_x = [];
            elements_y = [];
            elements_z = [];

            for m = 1:numSensors
                sensorCenter = [sensorLocations(m,1),sensorLocations(m,2),sensorLocations(m,3)];
                plot3(sensorCenter(1),sensorCenter(2),sensorCenter(3),'r.','MarkerSize',10);
                for k = 1:numElements

                    element_ids = fmdl.elems(k,:);

                    vertices = fmdl.nodes(element_ids,:);

                    [L,n] = computeIntegralAndNormalsInTetrahedron(vertices,sensorCenter);

                    S = zeros(3,1);
                    for j = 1:4
                        S = S + L(j)*n(:,j);
                    end

                    Rx(m,k) = dot(S,[1,0,0]);
                    Ry(m,k) = dot(S,[0,1,0]);
                    Rz(m,k) = dot(S,[0,0,1]);

                    % Stop whenever difference between analytical and
                    % numerical is too great ( problem must be in the
                    % analytical computaiton!, not the numerical!)

                    if abs(Rx(m,k)-RxNewtonCotes4(m,k))/abs(Rx(m,k))*100>1
                        patch('Vertices', vertices, 'Faces', [1 2 3;1 2 4;2 3 4;1 3 4], ...
                            'FaceColor', 'cyan', 'FaceAlpha', 0.5, 'EdgeColor', 'black');
                        axis equal
                        xlabel('X'); ylabel('Y'); zlabel('Z');
                        view(3)
                        grid on

                        countx = countx+1;
                        elements_x  = unique([elements_x,k]);
                    end

                    if abs(Ry(m,k)-RyNewtonCotes4(m,k))/abs(Ry(m,k))*100>1
                        patch('Vertices', vertices, 'Faces', [1 2 3;1 2 4;2 3 4;1 3 4], ...
                            'FaceColor', 'yellow', 'FaceAlpha', 0.5, 'EdgeColor', 'black');
                        axis equal
                        xlabel('X'); ylabel('Y'); zlabel('Z');
                        view(3)
                        grid on

                        county = county+1;
                        elements_y  = unique([elements_y,k]);
                    end

                    if abs(Rz(m,k)-RzNewtonCotes4(m,k))/abs(Rz(m,k))*100>1
                        patch('Vertices', vertices, 'Faces', [1 2 3;1 2 4;2 3 4;1 3 4], ...
                            'FaceColor', 'magenta', 'FaceAlpha', 0.5, 'EdgeColor', 'black');
                        axis equal
                        xlabel('X'); ylabel('Y'); zlabel('Z');
                        view(3)
                        grid on

                        countz = countz+1;
                        elements_z  = unique([elements_z,k]);
                    end

                end

            end


            errorsNewtonCotes4Rx = abs(Rx-RxNewtonCotes4)./abs(Rx)*100;
            errorsNewtonCotes4Rx = errorsNewtonCotes4Rx(:);

            errorsNewtonCotes4Ry = abs(Ry-RyNewtonCotes4)./abs(Ry)*100;
            errorsNewtonCotes4Ry = errorsNewtonCotes4Ry(:);

            errorsNewtonCotes4Rz = abs(Rz-RzNewtonCotes4)./abs(Rz)*100;
            errorsNewtonCotes4Rz = errorsNewtonCotes4Rz(:);

            % figure;
            %
            % hold on
            % plot(errorsNewtonCotes4Rx)
            % plot(errorsNewtonCotes4Ry)
            % plot(errorsNewtonCotes4Rz)
            % hold off;
            %
            % legend('x','y','z')
            %
            % ylim([0,1]);

            %% Test

            passedRx = all(errorsNewtonCotes4Rx<errorThresh);
            passedRy = all(errorsNewtonCotes4Ry<errorThresh);
            passedRz = all(errorsNewtonCotes4Rz<errorThresh);

            testCase.verifyTrue(passedRx,sprintf('Analytical and numerical Rx difference for N-C35 quadrature is less than %d (%%)\n',errorThresh));
            testCase.verifyTrue(passedRy,sprintf('Analytical and numerical Ry difference for N-C35 quadrature is less than %d (%%)\n',errorThresh));
            testCase.verifyTrue(passedRz,sprintf('Analytical and numerical Rz difference for N-C35 quadrature is less than %d (%%)\n',errorThresh));

            %% Print Results

            fprintf('%s\n', repmat('-', 1, 30));
            fprintf('%-10s\n', 'Mean, s.d. and max of relative errors (in %) for all entries of the Rx matrices:');

            fprintf('%-15s\n','Rx:');
            fprintf('%-15s%-15s%-25s%-10s\n','Quadrature', 'Mean (%)', 'Standard Deviation (%)','Maximum (%)');
            fprintf('%-15s%-15.9g%-25.9g%-10.9g\n','NewtonCotes35',mean(errorsNewtonCotes4Rx),std(errorsNewtonCotes4Rx),max(errorsNewtonCotes4Rx));

            fprintf('%s\n', repmat('-', 1, 30));

            fprintf('%-15s\n','Ry:');
            fprintf('%-15s%-15s%-25s%-10s\n','Quadrature', 'Mean (%)', 'Standard Deviation (%)','Maximum (%)');
            fprintf('%-15s%-15.9g%-25.9g%-10.9g\n','NewtonCotes35',mean(errorsNewtonCotes4Ry),std(errorsNewtonCotes4Ry),max(errorsNewtonCotes4Ry));
            fprintf('%s\n', repmat('-', 1, 30));

            fprintf('%-15s\n','Rz:');
            fprintf('%-15s%-15s%-25s%-10s\n','Quadrature', 'Mean (%)', 'Standard Deviation (%)','Maximum (%)');
            fprintf('%-15s%-15.9g%-25.9g%-10.9g\n','NewtonCotes35',mean(errorsNewtonCotes4Rz),std(errorsNewtonCotes4Rz),max(errorsNewtonCotes4Rz));
            fprintf('%s\n', repmat('-', 1, 30));

        end

        function testRMatrixDerivativeWithRespectToSensorPosition(testCase)
            %% Prepare workspace
            % Get the full path of the current script
            fullpath = mfilename('fullpath');
            % Extract just the folder
            script_folder = fileparts(fullpath);
            cd(script_folder);

            % Have to add the functions path manually so prepare_workspace runs
            parent_folder =fileparts(script_folder);
            addpath(genpath(fullfile(parent_folder,'functions')));

            model_folder = prepare_workspace(script_folder);

            rng(1);

            %% Model parameters
            z0 = 0.0058; %(Ohm m^2) is the contact impedance from the CEM article 58 Ohm cm^2
            l0 = 40e-3; %(m) the tank radius
            I0 = 2.4e-3;%(A) the magnitude of the injected current

            % The derived characteristic units
            V0 = z0*I0/(l0^2); %(V)
            sigma0 = l0/z0; %(S/m)
            J0 = I0/(l0^2);

            model_parameters = create_default_3d_model_parameters(l0, z0, sigma0, I0);

            model_parameters.maxsz = max(model_parameters.height,model_parameters.radius)/5;
            model_parameters.numOfElectrodesPerRing = 4;
            model_parameters.numOfSensors = 4;
            model_parameters.sensorRadius = model_parameters.radius*2;

            %% Make forward model

            [~,fmdls] = mk_mdeit_model(model_parameters,model_folder);

            fmdl = fmdls{1};

            sensorLocations = zeros(numel(fmdl.sensors),3);
            for i = 1: numel(fmdl.sensors)
                sensorLocations(i,:) = fmdl.sensors(i).position;
            end

            fprintf("Test Parameters: \n"+ ...
                "Error threshold: %g (%%)\n",testCase.testParameters.errorThresh);

            %% Compute Rx derivative with analytic expression integrated numerically

            dRpdx = dRmkj_dp(fmdl, 1, 1);
            dRpdy = dRmkj_dp(fmdl, 1, 2);
            dRpdz = dRmkj_dp(fmdl, 1, 3);
            %% Compute Rx derivative with respect to x-coordinate of sensor position with finite-differences method

            delta = 1e-9;
            [Rx_unperturbed,Ry_unperturbed,Rz_unperturbed] = ...
                testRMatrices.computeRmatricesNewtonCotes4(fmdl,sensorLocations);

            dRxdx_num = zeros(numel(fmdl.sensors),size(fmdl.elems,1));
            dRxdy_num = zeros(numel(fmdl.sensors),size(fmdl.elems,1));
            dRxdz_num = zeros(numel(fmdl.sensors),size(fmdl.elems,1));
            for l = 1:numel(fmdl.sensors)

                sensor_locations_perturbed = sensorLocations;
                sensor_locations_perturbed(l,1) = sensor_locations_perturbed(l,1) + delta;

                [Rx_perturbed,Ry_perturbed,Rz_perturbed] = ...
                    testRMatrices.computeRmatricesNewtonCotes4(fmdl,sensor_locations_perturbed);

                temp_x = (Rx_perturbed-Rx_unperturbed)/delta;
                temp_y = (Ry_perturbed-Ry_unperturbed)/delta;
                temp_z = (Rz_perturbed-Rz_unperturbed)/delta;

                dRxdx_num(l,:) = temp_x(l,:);
                dRxdy_num(l,:) = temp_y(l,:);
                dRxdz_num(l,:) = temp_z(l,:);
            end

            %% Compute Rx derivative with respect to x-coordinate of sensor position with complex step method


            function [dRpdx,dRpdy,dRpdz] = compute_dR_ds_complex(fmdl,sensor_locations_0,delta,p)
                n_stim = numel(fmdl.stimulation);
                n_sensors = numel(fmdl.sensors);
                n_elem = size(fmdl.elems,1);

                dRpdx = zeros(numel(fmdl.sensors),size(fmdl.elems,1));
                dRpdy = zeros(numel(fmdl.sensors),size(fmdl.elems,1));
                dRpdz = zeros(numel(fmdl.sensors),size(fmdl.elems,1));

                for m = 1:numel(fmdl.sensors)
                    fprintf('Done iteration %i of %i\n',m,numel(fmdl.sensors))

                    sensor_locations_perturbed = sensor_locations_0;
                    sensor_locations_perturbed(m,p) = sensor_locations_0(m,p) + 1i*delta;

                    [Rx_perturbed,Ry_perturbed,Rz_perturbed,~] = ...
                        compute_r_matrices_local(fmdl,sensor_locations_perturbed);

                    dRpdx(m,:) = imag(Rx_perturbed(m,:)) / delta;
                    dRpdy(m,:) =  imag(Ry_perturbed(m,:)) / delta;
                    dRpdz(m,:) = imag(Rz_perturbed(m,:)) / delta;
                end
            end

            p = 1;
            [dRpdx_c,dRpdy_c,dRpdz_c] = compute_dR_ds_complex(fmdl,sensorLocations,delta,p);


            %% Compare
            errorThresh = testCase.testParameters.errorThresh;

            % Test against finite differences
            error_x = (dRpdx(:)-dRxdx_num(:) )./ dRpdx(:)*100;
            error_y = (dRpdy(:)-dRxdy_num(:) )./ dRpdy(:)*100;
            error_z = (dRpdz(:)-dRxdz_num(:) )./ dRpdz(:)*100;

            % Test against complex step method
            error_x_c = (dRpdx(:)-dRpdx_c(:) )./ dRpdx(:)*100;
            error_y_c = (dRpdy(:)-dRpdy_c(:) )./ dRpdy(:)*100;
            error_z_c = (dRpdz(:)-dRpdz_c(:) )./ dRpdz(:)*100;

            testCase.verifyTrue(all(error_x<errorThresh),sprintf('Analytical and numerical dRxdx difference is more than %d (%%)\n',errorThresh));
            testCase.verifyTrue(all(error_y<errorThresh),sprintf('Analytical and numerical dRxdx difference is more than %d (%%)\n',errorThresh));
            testCase.verifyTrue(all(error_z<errorThresh),sprintf('Analytical and numerical dRxdx difference is more than %d (%%)\n',errorThresh));

            testCase.verifyTrue(all(error_x_c<errorThresh),sprintf('Analytical and numerical dRxdx difference is more than %d (%%)\n',errorThresh));
            testCase.verifyTrue(all(error_y_c<errorThresh),sprintf('Analytical and numerical dRxdx difference is more than %d (%%)\n',errorThresh));
            testCase.verifyTrue(all(error_z_c<errorThresh),sprintf('Analytical and numerical dRxdx difference is more than %d (%%)\n',errorThresh));

        end

        function testLambdaDerivativeWithRespectToSensorPosition(testCase)
            %% Prepare workspace
            % Get the full path of the current script
            fullpath = mfilename('fullpath');
            % Extract just the folder
            script_folder = fileparts(fullpath);
            cd(script_folder);

            % Have to add the functions path manually so prepare_workspace runs
            parent_folder =fileparts(script_folder);
            addpath(genpath(fullfile(parent_folder,'functions')));

            model_folder = prepare_workspace(script_folder);

            rng(1);

            %% Model parameters
            z0 = 0.0058; %(Ohm m^2) is the contact impedance from the CEM article 58 Ohm cm^2
            l0 = 40e-3; %(m) the tank radius
            I0 = 2.4e-3;%(A) the magnitude of the injected current

            % The derived characteristic units
            V0 = z0*I0/(l0^2); %(V)
            sigma0 = l0/z0; %(S/m)
            J0 = I0/(l0^2);

            model_parameters = create_default_3d_model_parameters(l0, z0, sigma0, I0);

            model_parameters.maxsz = max(model_parameters.height,model_parameters.radius)/5;
            model_parameters.numOfElectrodesPerRing = 4;
            model_parameters.numOfSensors = 4;
            model_parameters.sensorRadius = model_parameters.radius*2;
            model_parameters.material = struct();

            %% Make forward model

            [~,fmdls] = mk_mdeit_model(model_parameters,model_folder);

            fmdl = fmdls{1};

            sensorLocations = zeros(numel(fmdl.sensors),3);
            for i = 1: numel(fmdl.sensors)
                sensorLocations(i,:) = fmdl.sensors(i).position;
            end

            fprintf("Test Parameters: \n"+ ...
                "Error threshold: %g (%%)\n",testCase.testParameters.errorThresh);

            img = mk_image_mdeit(fmdl,sigma0);

            num_sensors = numel(fmdl.sensors);
            n_nodes = size(fmdl.nodes,1);
            %% Analitically with derived expression

            dlambda_analytical = compute_dlambda_ds(img,1,1);

            %% Numerically with finite differences
            delta = 1e-7;

            %function: Assign sensor locations
            function img = assign_sensor_locations(img,sensor_locations)
                assert(numel(img.fwd_model.sensors) == size(sensor_locations,1));

                for m = 1: numel(img.fwd_model.sensors)
                    img.fwd_model.sensors(m).position = sensor_locations(m,:);
                end
            end

            dlambda_numerical = zeros(n_nodes,num_sensors);

            lambda_unperturbed = compute_lambda_x(img);

            for l = 1:numel(img.fwd_model.sensors)

                sensor_locations_perturbed = sensorLocations;
                sensor_locations_perturbed(l,1) = sensor_locations_perturbed(l,1) + delta;

                img_perturbed = assign_sensor_locations(img,sensor_locations_perturbed);

                lambda_perturbed = compute_lambda_x(img_perturbed);

                dlambda_numerical(:,l) = (lambda_perturbed(:,l)-lambda_unperturbed(:,l))/delta;
            end

            %% Numerically with complex step method
            dlambda_complex = zeros(n_nodes,num_sensors);

            p=1;
            sensor_locations_0 = sensorLocations;
            for l = 1:numel(img.fwd_model.sensors)

                sensor_locations_perturbed = sensor_locations_0;
                sensor_locations_perturbed(l,p) = sensor_locations_0(l,p) + 1i*delta;
                img_perturbed = assign_sensor_locations(img,sensor_locations_perturbed);


                lambda_perturbed = compute_lambda_x(img_perturbed);

                dlambda_complex(:,l) = imag(lambda_perturbed(:,l)) / delta;
            end

            %% Compare
            errorThresh = testCase.testParameters.errorThresh;

            error_x_numerical = abs(dlambda_numerical(:)-dlambda_analytical(:))./ abs(dlambda_analytical(:))*100;
            error_x_complex = abs(dlambda_complex(:)-dlambda_analytical(:) )./ abs(dlambda_analytical(:))*100;

            testCase.verifyTrue(all(error_x_numerical<errorThresh),sprintf('Analytical and numerical difference is more than %d (%%)\n',errorThresh));
            testCase.verifyTrue(all(error_x_complex<errorThresh),sprintf('Analytical and numerical difference is more than %d (%%)\n',errorThresh));

        end

        function testGammaDerivativeWithRespectToSensorPosition(testCase)
            %% Prepare workspace
            % Get the full path of the current script
            fullpath = mfilename('fullpath');
            % Extract just the folder
            script_folder = fileparts(fullpath);
            cd(script_folder);

            % Have to add the functions path manually so prepare_workspace runs
            parent_folder =fileparts(script_folder);
            addpath(genpath(fullfile(parent_folder,'functions')));

            model_folder = prepare_workspace(script_folder);

            rng(1);

            %% Model parameters
            z0 = 0.0058; %(Ohm m^2) is the contact impedance from the CEM article 58 Ohm cm^2
            l0 = 40e-3; %(m) the tank radius
            I0 = 2.4e-3;%(A) the magnitude of the injected current

            % The derived characteristic units
            V0 = z0*I0/(l0^2); %(V)
            sigma0 = l0/z0; %(S/m)
            J0 = I0/(l0^2);

            model_parameters = create_default_3d_model_parameters(l0, z0, sigma0, I0);

            model_parameters.maxsz = max(model_parameters.height,model_parameters.radius)/5;
            model_parameters.numOfElectrodesPerRing = 4;
            model_parameters.numOfSensors = 4;
            model_parameters.sensorRadius = model_parameters.radius*2;
            model_parameters.material = struct();

            %% Make forward model

            [~,fmdls] = mk_mdeit_model(model_parameters,model_folder);

            fmdl = fmdls{1};

            sensorLocations = zeros(numel(fmdl.sensors),3);
            for i = 1: numel(fmdl.sensors)
                sensorLocations(i,:) = fmdl.sensors(i).position;
            end

            fprintf("Test Parameters: \n"+ ...
                "Error threshold: %g (%%)\n",testCase.testParameters.errorThresh);

            img = mk_image_mdeit(fmdl,sigma0);

            num_sensors = numel(fmdl.sensors);
            n_elems = size(fmdl.elems,1);
            n_nodes = size(fmdl.nodes,1);

            p=1;

            %% Numerically with finite differences
            delta = 1e-7;

            dgamma_numerical = zeros(num_sensors,n_nodes);

            img = compute_gamma_matrices_local(img);

            gamma_x_unperturbed = img.Gamma1;

            for m = 1:numel(img.fwd_model.sensors)

                sensor_locations_perturbed = sensorLocations;
                sensor_locations_perturbed(m,1) = sensor_locations_perturbed(m,1) + delta;

                img_perturbed = assign_sensor_locations(img,sensor_locations_perturbed);

                img_perturbed =  compute_gamma_matrices_local(img_perturbed);
                gamma_x_perturbed = img_perturbed.Gamma1;

                dgamma_numerical(m,:) = (gamma_x_perturbed(m,:)-gamma_x_unperturbed(m,:))/delta;
            end

            %% Numerically with complex step method
            dgamma_complex = zeros(num_sensors,n_nodes);

            for m = 1:numel(img.fwd_model.sensors)
                sensor_locations_perturbed = sensorLocations;
                sensor_locations_perturbed(m,p) = sensor_locations_perturbed(m,p) + 1i*delta;
                img_perturbed = assign_sensor_locations(img,sensor_locations_perturbed);

                img_perturbed =  compute_gamma_matrices_local(img_perturbed);
                gamma_x_perturbed = img_perturbed.Gamma1;

                dgamma_complex(m,:) = imag(gamma_x_perturbed(m,:)) / delta;
            end

            %% Compare
            errorThresh = testCase.testParameters.errorThresh;

            error_x = abs(dgamma_numerical(:)-dgamma_complex(:))./ abs(dgamma_complex(:))*100;

            testCase.verifyTrue(all(error_x<errorThresh),sprintf('Analytical and numerical difference is more than %d (%%)\n',errorThresh));

        end

        function testBDerivativeWithRespectToSensorPosition(testCase)
            %% Prepare workspace
            % Get the full path of the current script
            fullpath = mfilename('fullpath');
            % Extract just the folder
            script_folder = fileparts(fullpath);
            cd(script_folder);

            % Have to add the functions path manually so prepare_workspace runs
            parent_folder =fileparts(script_folder);
            addpath(genpath(fullfile(parent_folder,'functions')));

            model_folder = prepare_workspace(script_folder);

            rng(1);

            %% Model parameters
            z0 = 0.0058; %(Ohm m^2) is the contact impedance from the CEM article 58 Ohm cm^2
            l0 = 40e-3; %(m) the tank radius
            I0 = 2.4e-3;%(A) the magnitude of the injected current

            % The derived characteristic units
            V0 = z0*I0/(l0^2); %(V)
            sigma0 = l0/z0; %(S/m)
            J0 = I0/(l0^2);

            model_parameters = create_default_3d_model_parameters(l0, z0, sigma0, I0);

            model_parameters.maxsz = max(model_parameters.height,model_parameters.radius)/5;
            model_parameters.numOfElectrodesPerRing = 4;
            model_parameters.numOfSensors = 6;
            model_parameters.sensorRadius = model_parameters.radius*1.05;
            model_parameters.material = struct();

            %% Make forward model

            [~,fmdls] = mk_mdeit_model(model_parameters,model_folder);

            fmdl = fmdls{1};

            sensorLocations = zeros(numel(fmdl.sensors),3);
            for i = 1: numel(fmdl.sensors)
                sensorLocations(i,:) = fmdl.sensors(i).position;
            end

            fprintf("Test Parameters: \n"+ ...
                "Error threshold: %g (%%)\n",testCase.testParameters.errorThresh);

            img = mk_image_mdeit(fmdl,sigma0);


            n_sensors = numel(img.fwd_model.sensors);
            n_stim = numel(img.fwd_model.stimulation);
            n_elem = size(img.fwd_model.elems,1);
            n_nodes = size(fmdl.nodes,1);
            %%

            function img = assign_sensor_locations(img,sensor_locations)
                assert(numel(img.fwd_model.sensors) == size(sensor_locations,1));

                for l = 1: numel(img.fwd_model.sensors)
                    img.fwd_model.sensors(l).position = sensor_locations(l,:);
                end
            end

            %% Compute dB analitically

            dim = 1; %x
            p=1 ; % wrt sensor x-coordinate

            dB_analytical = compute_dB_ds(img,dim,p);

            %% Compute dB numerically

            delta = 1e-7;

            % Compute EIT forward solution for each current injection pattern
            u = fwd_solve(img);
            u = u.volt;

            img = compute_gamma_matrices_local(img);

            dB_numerical = zeros(numel(img.fwd_model.sensors),n_stim);

            Bx_unperturbed = img.Gamma1*u;
            for m = 1:numel(img.fwd_model.sensors)
                fprintf('Done iteration %i of %i\n',m,numel(img.fwd_model.sensors))

                sensor_locations_perturbed = sensorLocations;
                sensor_locations_perturbed(m,p) = sensor_locations_perturbed(m,p) + delta;
                img_perturbed = assign_sensor_locations(img,sensor_locations_perturbed);

                img_perturbed = compute_gamma_matrices_local(img_perturbed);

                Bx_perturbed = img_perturbed.Gamma1*u;

                dB_numerical(m,:) = (Bx_perturbed(m,:)-Bx_unperturbed(m,:)) / delta;
            end

            %% Compare
            errorThresh = testCase.testParameters.errorThresh;

            error = abs(dB_analytical(:)-dB_numerical(:))./abs(dB_analytical(:))*100;

            testCase.verifyTrue(all(error<errorThresh),sprintf('Analytical and numerical difference is more than %d (%%)\n',errorThresh));

        end

        function testJacobianDerivativeWithRespectToSensorPosition(testCase)
            %% Prepare workspace
            % Get the full path of the current script
            fullpath = mfilename('fullpath');
            % Extract just the folder
            script_folder = fileparts(fullpath);
            cd(script_folder);

            % Have to add the functions path manually so prepare_workspace runs
            parent_folder =fileparts(script_folder);
            addpath(genpath(fullfile(parent_folder,'functions')));

            model_folder = prepare_workspace(script_folder);

            rng(1);

            %% Model parameters
            z0 = 0.0058; %(Ohm m^2) is the contact impedance from the CEM article 58 Ohm cm^2
            l0 = 40e-3; %(m) the tank radius
            I0 = 2.4e-3;%(A) the magnitude of the injected current

            % The derived characteristic units
            V0 = z0*I0/(l0^2); %(V)
            sigma0 = l0/z0; %(S/m)
            J0 = I0/(l0^2);

            model_parameters = create_default_3d_model_parameters(l0, z0, sigma0, I0);

            model_parameters.maxsz = max(model_parameters.height,model_parameters.radius)/5;
            model_parameters.numOfElectrodesPerRing = 4;
            model_parameters.numOfSensors = 6;
            model_parameters.sensorRadius = model_parameters.radius*1.5;
            model_parameters.material = struct();

            %% Make forward model

            [~,fmdls] = mk_mdeit_model(model_parameters,model_folder);

            fmdl = fmdls{1};

            sensorLocations = zeros(numel(fmdl.sensors),3);
            for i = 1: numel(fmdl.sensors)
                sensorLocations(i,:) = fmdl.sensors(i).position;
            end

            fprintf("Test Parameters: \n"+ ...
                "Error threshold: %g (%%)\n",testCase.testParameters.errorThresh);

            img = mk_image_mdeit(fmdl,sigma0);

            n_sensors = numel(img.fwd_model.sensors);
            n_stim = numel(img.fwd_model.stimulation);
            n_elem = size(img.fwd_model.elems,1);
            n_nodes = size(fmdl.nodes,1);


            %% Compute ....


            % BUG!!!!!!!!!! For some reason, can't figure out why, finite
            % differences is not giving the same result as complex step and
            % analytical derivative:
            % Solution: Took me a while, but the problem was that sometimes
            % I was taking the adjoint transposed with ('), instead of just
            % the transposed (.'). That made some things change sign with
            % the complex step method! I noticed this because the FD method
            % was equal to the analytic when dJ = dJ1-dJ2 and the complex
            % step method was equal to the analytic when dJ = dJ1 + dJ2;


            A = @(x) M(img,x);

            delta = 1e-3;

            for dim = 1:3
                for p = 1:3
                    % Analitically with expression
                    [dJ_analytical,~,~] = compute_dJ_ds(img,dim,p);

                    % Numerically with complex step method
                    dJ_complex =  compute_dJ_ds_complex(img,A,sensorLocations,delta,dim,p);

                    % Numerically with central differences
                    dJ_numerical =  compute_dJ_ds_central_differences(img,A,sensorLocations,delta,dim,p);

                    % Compare

                    err_complex_fd = norm(dJ_complex(:) - dJ_numerical(:)) / norm(dJ_complex(:));
                    err_analytical_fd = norm(dJ_analytical(:) - dJ_numerical(:)) / norm(dJ_analytical(:));
                    err_analyctical_complex = norm(dJ_analytical(:) - dJ_complex(:)) / norm(dJ_complex(:));

                    errorThresh = testCase.testParameters.errorThresh;

                    testCase.verifyTrue(all(err_complex_fd<errorThresh),sprintf('Analytical and numerical difference is more than %d (%%)\n',errorThresh));
                    testCase.verifyTrue(all(err_analytical_fd<errorThresh),sprintf('Analytical and numerical difference is more than %d (%%)\n',errorThresh));
                    testCase.verifyTrue(all(err_analyctical_complex<errorThresh),sprintf('Analytical and numerical difference is more than %d (%%)\n',errorThresh));

                end
            end

            % finish
        end
       
        function testCostFunctionGradientWithRespectToSensorPosition(testCase)
            %% Prepare workspace
            % Get the full path of the current script
            fullpath = mfilename('fullpath');
            % Extract just the folder
            script_folder = fileparts(fullpath);
            cd(script_folder);

            % Have to add the functions path manually so prepare_workspace runs
            parent_folder =fileparts(script_folder);
            addpath(genpath(fullfile(parent_folder,'functions')));

            model_folder = prepare_workspace(script_folder);

            rng(1);

            %% Model parameters
            z0 = 0.0058; %(Ohm m^2) is the contact impedance from the CEM article 58 Ohm cm^2
            l0 = 40e-3; %(m) the tank radius
            I0 = 2.4e-3;%(A) the magnitude of the injected current

            % The derived characteristic units
            V0 = z0*I0/(l0^2); %(V)
            sigma0 = l0/z0; %(S/m)
            J0 = I0/(l0^2);

            model_parameters = create_default_3d_model_parameters(l0, z0, sigma0, I0);

            model_parameters.maxsz = max(model_parameters.height,model_parameters.radius)/5;
            model_parameters.numOfElectrodesPerRing = 4;
            model_parameters.numOfSensors = 2;
            model_parameters.sensorRadius = model_parameters.radius*1.5;
            model_parameters.material.type = 'spherical';

            background_conductivity = 3.28e-1/sigma0;
            anomaly_conductivity = background_conductivity/10;
            %% Make forward model

            [~,fmdls] = mk_mdeit_model(model_parameters,model_folder);

            fmdl = fmdls{1};

            sensorLocations = zeros(numel(fmdl.sensors),3);
            for i = 1: numel(fmdl.sensors)
                sensorLocations(i,:) = fmdl.sensors(i).position;
            end

            fprintf("Test Parameters: \n"+ ...
                "Error threshold: %g (%%)\n",testCase.testParameters.errorThresh);

            % Make homogeneous image
            imgh = mk_image_mdeit(fmdl,background_conductivity);

            % Add plastic cylinder
            img = add_material_properties(imgh, [background_conductivity,anomaly_conductivity]);

            show_fem(img);

            n_sensors = numel(img.fwd_model.sensors);
            n_stim = numel(img.fwd_model.stimulation);
            n_elem = size(img.fwd_model.elems,1);
            n_nodes = size(fmdl.nodes,1);

            %% Assume white noise and white prior

            B = fwd_solve_mdeit(img);
            max_B = max([abs(B.Bx(:));abs(B.By(:));abs(B.Bz(:))]);

            noise_std_deviation = max_B/20;
            prior_std_deviation = max(img.elem_data)/10;

            % White noise with zero mean and \mu variance
            Gamma_noise = noise_std_deviation^2.*speye(n_stim*n_sensors,n_stim*n_sensors);

            % Gaussian prior with mean = mean of initial conductivity and
            % prior_std_deviation

            % prior_mean = mean(img.elem_data);
            Gamma_prior = prior_std_deviation^2.*speye(n_elem,n_elem);

            %% Compare analytical and numerical computation of gradient
            A = @(x) M(img,x);
            delta = 1e-5;

            errorThresh = testCase.testParameters.errorThresh;

            function dphidp_numerical = compute_cost_function_gradient_central_differences(img,A,dim,p,Gamma_prior,Gamma_noise,delta)

                n_sensors = numel(img.fwd_model.sensors);
                dphidp_numerical = zeros(1,n_sensors);

                sensor_locations_0_local = fetch_sensor_locations(img);

                for sensor_id = 1:n_sensors
                    sensor_locations_local = sensor_locations_0_local;
                    sensor_locations_local(sensor_id,p) = sensor_locations_local(sensor_id,p) + delta;

                    % assign sensor locations is done inside
                    % compute_cost_function
                    cost_m_plus = compute_cost_function(img,sensor_locations_local,Gamma_prior,Gamma_noise,A,dim);

                    sensor_locations_local = sensor_locations_0_local;
                    sensor_locations_local(sensor_id,p) = sensor_locations_local(sensor_id,p) - delta;

                    cost_m_minus = compute_cost_function(img,sensor_locations_local,Gamma_prior,Gamma_noise,A,dim);

                    dphidp_numerical(sensor_id) = (cost_m_plus-cost_m_minus)/(2*delta);
                end
            end

            for dim = 1:3
                for p = 1:3

                    %Compute the gradient of the cost function analytically
                    dphidp = compute_cost_function_gradient(...
                        img,sensorLocations,Gamma_prior,Gamma_noise,A,dim,p);

                    %Compute the gradient of the cost function numerically with central-differences
                    dphidp_numerical = compute_cost_function_gradient_central_differences(img,A,dim,p,Gamma_prior,Gamma_noise,delta);

                    err = norm(dphidp_numerical(:) - dphidp(:)) / norm(dphidp(:));
                    testCase.verifyTrue(all(err<errorThresh),sprintf('Analytical and numerical difference is more than %d (%%)\n',errorThresh));
                end
            end

        end

        function testAOptimality(testCase)
            
            %% Prepare workspace
            % Get the full path of the current script
            fullpath = mfilename('fullpath');
            % Extract just the folder
            script_folder = fileparts(fullpath);
            cd(script_folder);

            % Have to add the functions path manually so prepare_workspace runs
            parent_folder =fileparts(script_folder);
            addpath(genpath(fullfile(parent_folder,'functions')));

            model_folder = prepare_workspace(script_folder);

            rng(1);

            %% Model parameters
            z0 = 0.0058; %(Ohm m^2) is the contact impedance from the CEM article 58 Ohm cm^2
            l0 = 40e-3; %(m) the tank radius
            I0 = 2.4e-3;%(A) the magnitude of the injected current

            % The derived characteristic units
            V0 = z0*I0/(l0^2); %(V)
            sigma0 = l0/z0; %(S/m)
            J0 = I0/(l0^2);

            model_parameters = create_default_3d_model_parameters(l0, z0, sigma0, I0);

            model_parameters.maxsz = max(model_parameters.height,model_parameters.radius)/5;
            model_parameters.numOfElectrodesPerRing = 4;
            model_parameters.numOfSensors = 4;
            model_parameters.sensorRadius = model_parameters.radius*1.5;
            model_parameters.material.type = 'spherical';

            background_conductivity = 3.28e-1/sigma0;
            anomaly_conductivity = background_conductivity/10;
            %% Make forward model

            [~,fmdls] = mk_mdeit_model(model_parameters,model_folder);

            fmdl = fmdls{1};

            sensorLocations = zeros(numel(fmdl.sensors),3);
            for i = 1: numel(fmdl.sensors)
                sensorLocations(i,:) = fmdl.sensors(i).position;
            end

            fprintf("Test Parameters: \n"+ ...
                "Error threshold: %g (%%)\n",testCase.testParameters.errorThresh);

            % Make homogeneous image
            imgh = mk_image_mdeit(fmdl,background_conductivity);

            % Add plastic cylinder
            img = add_material_properties(imgh, [background_conductivity,anomaly_conductivity]);

            show_fem(img);

            n_sensors = numel(img.fwd_model.sensors);
            n_stim = numel(img.fwd_model.stimulation);
            n_elem = size(img.fwd_model.elems,1);
            n_nodes = size(fmdl.nodes,1);

            %% Assume white noise and white prior

            B = fwd_solve_mdeit(img);
            max_B = max([abs(B.Bx(:));abs(B.By(:));abs(B.Bz(:))]);

            noise_std_deviation = max_B/20;
            prior_std_deviation = max(img.elem_data)/10;

            % White noise with zero mean and \mu variance
            Gamma_noise = noise_std_deviation^2.*speye(n_stim*n_sensors,n_stim*n_sensors);

            % Gaussian prior with mean = mean of initial conductivity and
            % prior_std_deviation

            % prior_mean = mean(img.elem_data);
            Gamma_prior = prior_std_deviation^2.*speye(n_elem,n_elem);
            
            %% Perform optimization with lbfgs
            
            A = @(x) M(img,x);
            dim = 1;

            % Gamma_prior is constant
            % Gamma_noise is constant
            % A is constant
            % dim is constant
            % img is constant, gets changed inside compute_cost_function
            function out = f(img,x,Gamma_prior,Gamma_noise,A,dim)

                sensor_locations = vector_to_sensor_locations(x);
               
                out = compute_cost_function(img,sensor_locations,Gamma_prior,Gamma_noise,A,dim);
            
            end
            
            function dphi = g(img,x,Gamma_prior,Gamma_noise,A,dim)

                n_sensors = numel(img.fwd_model.sensors);

                sensor_locations = vector_to_sensor_locations(x);
                
                dphi = zeros(n_sensors*3,1);
                for p = 1:3
                    dphidp = compute_cost_function_gradient(...
                        img,sensor_locations,Gamma_prior,Gamma_noise,A,dim,p);
                    ids = (p-1)*n_sensors+1:p*n_sensors;
                    dphi(ids) = dphidp;
                end

            end
            

            func = @(x) f(img,x,Gamma_prior,Gamma_noise,A,dim);
            grad = @(x) g(img,x,Gamma_prior,Gamma_noise,A,dim);
            
            % Test if the functions are working
            sensor_locations_0 = fetch_sensor_locations(img);
            x0 = sensor_locations_to_vector(sensor_locations_0);

            out_f = func(x0);           
            out_g = grad(x0);
            
            %% Optimize with L-BFGS
            tol = 1e-9;
            max_iterations = 10;
            options.display = true;
            options.maxIterations = 10;
            options.tol = tol;
            % [xk,k] = LBFGS(func,grad,x0,options);
            
            % [xk,k] = BFGS(func,grad,x0,1e-9,10)
        
            [xk,k] = steepestDescent(func,grad,x0,1e-9,max_iterations);


            sensor_locations_k = vector_to_sensor_locations(xk);
            
            figure
            show_fem(img);
            plot_sensors(img);
            img = assign_sensor_locations(img,sensor_locations_k);
            plot_sensors(img,false,'r');
            box on;grid on;
            drawnow
            

            % %% Compute and visualize the cost function in a ball around the initial p-coordinates of the 2 sensors
            % 
            % sensor_p_coordinates = zeros(n_sensors,1);
            % for m = 1:n_sensors
            %     sensor_p_coordinates(m) = img.fwd_model.sensors(m).position(1);
            % end
            % 
            % sensor_locations_0 = fetch_sensor_locations(img);
            % 
            % 
            % sensor_radius =  sqrt(sensor_locations_0(:,1).^2+sensor_locations_0(:,2).^2);
            % 
            % R = 0.05*(min(sensor_radius)-model_parameters.radius);
            % 
            % % Number of different positions for each sensor to plot the
            % % cost function at
            % N  = 5;
            % 
            % sensor_locations_around_ball = zeros(n_sensors,3,N);
            % 
            % for m = 1:n_sensors
            %     for n = 1:3
            %         if n==p
            %             sensor_locations_around_ball(m,n,:) = linspace(sensor_p_coordinates(m)-R,sensor_p_coordinates(m)+R,N);
            %         else
            %             sensor_locations_around_ball(m,n,:) = ones(N,1)*sensor_locations_0(m,p);
            %         end
            %     end
            % end
            % 
            % % Include the starting sensor position in the list of sensor
            % % positions
            % % sensor_locations_around_ball = cat(3,sensor_locations_around_ball,sensor_locations_0);
            % % sensor_locations_around_ball = sort(sensor_locations_around_ball,3);
            % 
            % [X1,X2] = meshgrid(sensor_locations_around_ball(1,p,:),sensor_locations_around_ball(2,p,:));
            % 
            % % Just a visualization
            % figure
            % hold on
            % show_fem(img);
            % for i = 1:size(sensor_locations_around_ball,3)
            %     sensor_locations = sensor_locations_0;
            %     sensor_locations(1,p) = X1(i,i);
            %     sensor_locations(2,p) = X2(i,i);
            %     img_temp = assign_sensor_locations(img,sensor_locations);
            %     plot_sensors(img_temp)
            % end
            % hold off
            % drawnow;
            % 
            % % Compute cost function at different p-coordinates
            % cost_box = zeros(N,N);
            % 
            % for i = 1:size(sensor_locations_around_ball,3)
            %     for j = 1:size(sensor_locations_around_ball,3)
            %         fprintf('Evaluating cost (%i/%i)\n',...
            %             size(sensor_locations_around_ball,3)*(i-1)+j,size(sensor_locations_around_ball,3)^2)
            % 
            %         sensor_locations = sensor_locations_0;
            % 
            %         sensor_locations(1,p) = X1(i,j);
            %         sensor_locations(2,p) = X2(i,j);
            % 
            %         cost_box(i,j) = compute_cost_function(img,sensor_locations,Gamma_prior,Gamma_noise,A,dim);
            %     end
            % end
            % 
            % cost_at_initial_position = compute_cost_function(img,sensor_locations_0,Gamma_prior,Gamma_noise,A,dim);
            % 
            % % Estimate gradient
            % % Compute spacing in each direction (assuming uniform grid)
            % dx = X1(1,2) - X1(1,1);  % spacing along X1 (columns)
            % dy = X2(2,1) - X2(1,1);  % spacing along X2 (rows)
            % 
            % % Compute numerical gradient
            % [phi_x, phi_y] = gradient(cost_box, dx, dy);
            % %% Plot cost function
            % norm_dphip = sqrt(dphidp(1)^2+dphidp(2)^2);
            % norm_phi = sqrt(phi_x.^2+phi_y.^2);
            % scale = abs(max(X1(:))-min(X1(:)));
            % figure
            % hold on
            % quiver3(sensor_locations_0(1,1),sensor_locations_0(2,1),cost_at_initial_position,...
            %     dphidp(1)/norm_dphip,dphidp(2)/norm_dphip,0,...
            %     scale,'LineStyle','-');
            % quiver3(X1,X2,cost_box,...
            %     phi_x./norm_phi ,phi_y./norm_phi ,zeros(size(phi_x)),...
            %     10*scale,'LineStyle','--')
            % plot3(sensor_locations_0(1,1),sensor_locations_0(2,1),cost_at_initial_position,'r.','MarkerSize',10)
            % mesh(X1,X2,cost_box);
            % xlabel('$p$-coordinate of sensor $1$','Interpreter','latex')
            % ylabel('$p$-coordinate of sensor $2$','Interpreter','latex')
            % box on;grid on;grid minor;
            % axis('square')
            % hold off
        end

    end


    %% Static methods
    methods (Static, Access = private)
        function [Rx, Ry, Rz] = computeRmatricesKeast0(fmdl, sensorLocations)
            numElements = size(fmdl.elems,1);
            numSensors  = size(sensorLocations,1);

            Rx = zeros(numSensors, numElements);
            Ry = zeros(numSensors, numElements);
            Rz = zeros(numSensors, numElements);

            % Keast order-0 quadrature point and weight
            coord = [0.25; 0.25; 0.25];
            weight = 1.0;

            for m = 1:numSensors
                rm = sensorLocations(m,:);

                % Define the function to integrate
                fun = @(x, y, z) (rm - [x, y, z]) ./ (sum((rm - [x, y, z]).^2)^1.5);

                for k = 1:numElements
                    v = fmdl.nodes(fmdl.elems(k,:),:);

                    % Jacobian matrix and determinant
                    J = [(v(2,:)-v(1,:))', (v(3,:)-v(1,:))', (v(4,:)-v(1,:))'];
                    detJ = abs(det(J));

                    % Quadrature point in physical tetrahedron
                    xi = v(1,:)' + J * coord;

                    % Evaluate integrand
                    R = weight * fun(xi(1), xi(2), xi(3));

                    % Scale by physical volume (1/6 from reference tetrahedron)
                    R = (detJ / 6) * R;

                    % Assign
                    Rx(m,k) = R(1);
                    Ry(m,k) = R(2);
                    Rz(m,k) = R(3);
                end
            end
        end

        function [Rx,Ry,Rz] = computeRmatricesNewtonCotes1(fmdl,sensorLocations)
            numElements = size(fmdl.elems,1);
            numSensors = size(sensorLocations,1);

            Rx = zeros(numSensors,numElements);
            Ry = zeros(numSensors,numElements);
            Rz = zeros(numSensors,numElements);

            coord = [0.0000000000000000  0.0000000000000000  0.0000000000000000
                0.0000000000000000  0.0000000000000000  1.0000000000000000
                0.0000000000000000  1.0000000000000000  0.0000000000000000
                1.0000000000000000  0.0000000000000000  0.0000000000000000];

            weights = [  0.2500000000000000
                0.2500000000000000
                0.2500000000000000
                0.2500000000000000];

            for m = 1:numSensors
                rm = sensorLocations(m,:);
                fun = @(x,y,z) (rm - [x,y,z])./ (sum((rm - [x, y, z]).^2)^1.5);
                for k = 1:numElements

                    %Find the vertices of the tetrahedron
                    v = fmdl.nodes(fmdl.elems(k,:),:);

                    J = [(v(2,:)-v(1,:))',(v(3,:)-v(1,:))',(v(4,:)-v(1,:))'];

                    detJ = abs(det(J));

                    R = 0;
                    for i = 1:length(weights)

                        r = coord(i,1);
                        s = coord(i,2);
                        t = coord(i,3);

                        xi = v(1,:)' + J * [r; s; t];

                        R = R + weights(i)*fun(xi(1),xi(2),xi(3));
                    end

                    R =   (detJ / 6) * R;

                    Rx(m,k) = R(1);
                    Ry(m,k) = R(2);
                    Rz(m,k) = R(3);
                end
            end
        end

        function [Rx,Ry,Rz] = computeRmatricesNewtonCotes2(fmdl,sensorLocations)
            numElements = size(fmdl.elems,1);
            numSensors = size(sensorLocations,1);

            Rx = zeros(numSensors,numElements);
            Ry = zeros(numSensors,numElements);
            Rz = zeros(numSensors,numElements);

            coord = [  0.0000000000000000  0.0000000000000000  0.0000000000000000
                0.0000000000000000  0.0000000000000000  1.0000000000000000
                0.0000000000000000  1.0000000000000000  0.0000000000000000
                1.0000000000000000  0.0000000000000000  0.0000000000000000
                0.5000000000000000  0.5000000000000000  0.0000000000000000
                0.5000000000000000  0.0000000000000000  0.5000000000000000
                0.5000000000000000  0.0000000000000000  0.0000000000000000
                0.0000000000000000  0.5000000000000000  0.5000000000000000
                0.0000000000000000  0.5000000000000000  0.0000000000000000
                0.0000000000000000  0.0000000000000000  0.5000000000000000];

            weights = [ -0.0500000000000000
                -0.0500000000000000
                -0.0500000000000000
                -0.0500000000000000
                0.2000000000000000
                0.2000000000000000
                0.2000000000000000
                0.2000000000000000
                0.2000000000000000
                0.2000000000000000];

            for m = 1:numSensors
                rm = sensorLocations(m,:);
                fun = @(x,y,z) (rm - [x,y,z])./ (sum((rm - [x, y, z]).^2)^1.5);
                for k = 1:numElements

                    %Find the vertices of the tetrahedron
                    v = fmdl.nodes(fmdl.elems(k,:),:);

                    J = [(v(2,:)-v(1,:))',(v(3,:)-v(1,:))',(v(4,:)-v(1,:))'];

                    detJ = abs(det(J));

                    R = 0;
                    for i = 1:length(weights)

                        r = coord(i,1);
                        s = coord(i,2);
                        t = coord(i,3);

                        xi = v(1,:)' + J * [r; s; t];

                        R = R + weights(i)*fun(xi(1),xi(2),xi(3));
                    end

                    R =  (detJ / 6) * R;
                    Rx(m,k) = R(1);
                    Ry(m,k) = R(2);
                    Rz(m,k) = R(3);
                end
            end
        end

        function [Rx,Ry,Rz] = computeRmatricesNewtonCotes4(fmdl,sensorLocations)
            numElements = size(fmdl.elems,1);
            numSensors = size(sensorLocations,1);

            Rx = zeros(numSensors,numElements);
            Ry = zeros(numSensors,numElements);
            Rz = zeros(numSensors,numElements);

            coord = [    0.0000000000000000  0.0000000000000000  0.0000000000000000
                0.0000000000000000  0.0000000000000000  1.0000000000000000
                0.0000000000000000  1.0000000000000000  0.0000000000000000
                1.0000000000000000  0.0000000000000000  0.0000000000000000
                0.0000000000000000  0.0000000000000000  0.7500000000000000
                0.0000000000000000  0.0000000000000000  0.2500000000000000
                0.0000000000000000  0.7500000000000000  0.0000000000000000
                0.0000000000000000  0.7500000000000000  0.2500000000000000
                0.0000000000000000  0.2500000000000000  0.0000000000000000
                0.0000000000000000  0.2500000000000000  0.7500000000000000
                0.7500000000000000  0.0000000000000000  0.0000000000000000
                0.7500000000000000  0.0000000000000000  0.2500000000000000
                0.7500000000000000  0.2500000000000000  0.0000000000000000
                0.2500000000000000  0.0000000000000000  0.0000000000000000
                0.2500000000000000  0.0000000000000000  0.7500000000000000
                0.2500000000000000  0.7500000000000000  0.0000000000000000
                0.5000000000000000  0.5000000000000000  0.0000000000000000
                0.5000000000000000  0.0000000000000000  0.5000000000000000
                0.5000000000000000  0.0000000000000000  0.0000000000000000
                0.0000000000000000  0.5000000000000000  0.5000000000000000
                0.0000000000000000  0.5000000000000000  0.0000000000000000
                0.0000000000000000  0.0000000000000000  0.5000000000000000
                0.2500000000000000  0.2500000000000000  0.0000000000000000
                0.2500000000000000  0.2500000000000000  0.5000000000000000
                0.2500000000000000  0.0000000000000000  0.2500000000000000
                0.2500000000000000  0.0000000000000000  0.5000000000000000
                0.2500000000000000  0.5000000000000000  0.2500000000000000
                0.2500000000000000  0.5000000000000000  0.0000000000000000
                0.0000000000000000  0.2500000000000000  0.2500000000000000
                0.0000000000000000  0.2500000000000000  0.5000000000000000
                0.0000000000000000  0.5000000000000000  0.2500000000000000
                0.5000000000000000  0.2500000000000000  0.2500000000000000
                0.5000000000000000  0.2500000000000000  0.0000000000000000
                0.5000000000000000  0.0000000000000000  0.2500000000000000
                0.2500000000000000  0.2500000000000000  0.2500000000000000];

            weights = [  -0.0119047619047619
                -0.0119047619047619
                -0.0119047619047619
                -0.0119047619047619
                0.0380952380952381
                0.0380952380952381
                0.0380952380952381
                0.0380952380952381
                0.0380952380952381
                0.0380952380952381
                0.0380952380952381
                0.0380952380952381
                0.0380952380952381
                0.0380952380952381
                0.0380952380952381
                0.0380952380952381
                -0.0285714285714286
                -0.0285714285714286
                -0.0285714285714286
                -0.0285714285714286
                -0.0285714285714286
                -0.0285714285714286
                0.0380952380952381
                0.0380952380952381
                0.0380952380952381
                0.0380952380952381
                0.0380952380952381
                0.0380952380952381
                0.0380952380952381
                0.0380952380952381
                0.0380952380952381
                0.0380952380952381
                0.0380952380952381
                0.0380952380952381
                0.3047619047619048];

            for m = 1:numSensors
                rm = sensorLocations(m,:);
                fun = @(x,y,z) (rm - [x,y,z])./ (sum((rm - [x, y, z]).^2)^1.5);
                for k = 1:numElements

                    %Find the vertices of the tetrahedron
                    v = fmdl.nodes(fmdl.elems(k,:),:);

                    J = [(v(2,:)-v(1,:))',(v(3,:)-v(1,:))',(v(4,:)-v(1,:))'];

                    detJ = abs(det(J));

                    R = 0;
                    for i = 1:length(weights)

                        r = coord(i,1);
                        s = coord(i,2);
                        t = coord(i,3);

                        xi = v(1,:)' + J * [r; s; t];

                        R = R + weights(i)*fun(xi(1),xi(2),xi(3));
                    end

                    R =  (detJ / 6) * R;
                    Rx(m,k) = R(1);
                    Ry(m,k) = R(2);
                    Rz(m,k) = R(3);
                end
            end
        end

    end


end

function out = M(img,sigma)

numNodes = size(img.fwd_model.nodes,1);

img.elem_data = sigma;
s_mat = system_mat_1st_order(img);

Ac = s_mat.E(1:numNodes,1:numNodes);
Ae = s_mat.E(1:numNodes,numNodes+1:end);
Ad = s_mat.E(numNodes+1:end,numNodes+1:end);

out = Ac-Ae*inv(Ad)*Ae';
end

function lambda = compute_lambda_x(img)

n_nodes = size(img.fwd_model.nodes,1);
num_sensors = numel(img.fwd_model.sensors);

sigma = img.elem_data;

% Compute Gamma matrices
img = compute_gamma_matrices_local(img);
Gamma = img.Gamma1;

% Solve the adjoint problem for each sensor to get lambda vectors
lambda = zeros(n_nodes,num_sensors);

A_matrix = M(img,sigma);

% Jacobi preconditioner - matrix free
d = sqrt(diag(A_matrix));        % vector of diagonal entries

Mfun = @(x) x ./ d;              % left preconditioner  M^{-1} x
Nfun = @(x) x ./ d;              % right preconditioner

GammaT = Gamma.';

parfor m = 1:num_sensors
    [lambda(:,m),~,~] = pcg(A_matrix,-GammaT(:,m),1e-6,numel(sigma),Mfun,Nfun);
end

end

function img = compute_gamma_matrices_local(img)

mu_factor = img.fwd_model.mu0/(4*pi);

num_sensors = numel(img.fwd_model.sensors);

sensor_locations = zeros(numel(img.fwd_model.sensors),3);

for i = 1: numel(img.fwd_model.sensors)
    sensor_locations(i,:) = img.fwd_model.sensors(i).position;
end

[Rx,Ry,Rz,fmdl] = compute_r_matrices_local(img.fwd_model,sensor_locations);
img.fwd_model = fmdl; % THIS IS CRUCIAL!!!!! so the R matrices are updated whenever compute_gamma_matrices_local is called, otherwise they will be the same as the initial img.fwd_model, which might be wrong!

% Convenience handles
R.Rx = Rx;
R.Ry = Ry;
R.Rz = Rz;
G = img.fwd_model.G;

% Sigma = sparse(1:length(img.elem_data), 1:length(img.elem_data), img.elem_data);
Sigma = spdiags(img.elem_data(:), 0, length(img.elem_data), length(img.elem_data));

% NEW: The matrix g_{dl}^m, gives the components of the measurement axis of
% sensor m on the canonical R^3 basis

g = zeros(num_sensors,3,3);
for m = 1:numel(img.fwd_model.sensors)
    g(m,:,:) = [...
        img.fwd_model.sensors(m).axes.axis1;
        img.fwd_model.sensors(m).axes.axis2;
        img.fwd_model.sensors(m).axes.axis3];
end

Cx = ( -R.Rz * Sigma * G.Gy +  R.Ry * Sigma * G.Gz );
Cy = ( -R.Rx * Sigma * G.Gz +  R.Rz * Sigma * G.Gx );
Cz = ( -R.Ry * Sigma * G.Gx +  R.Rx * Sigma * G.Gy );

Gamma1 = mu_factor*(g(:,1,1).*Cx + g(:,1,2).*Cy + g(:,1,3).*Cz);
Gamma2 = mu_factor*(g(:,2,1).*Cx + g(:,2,2).*Cy + g(:,2,3).*Cz);
Gamma3 = mu_factor*(g(:,3,1).*Cx + g(:,3,2).*Cy + g(:,3,3).*Cz);

img.Gamma1 = Gamma1;
img.Gamma2 = Gamma2;
img.Gamma3 = Gamma3;

end

function [Rx,Ry,Rz,fmdl] = compute_r_matrices_local(fmdl,sensor_locations)
numElements = size(fmdl.elems,1);
numSensors = size(sensor_locations,1);

Rx = zeros(numSensors,numElements);
Ry = zeros(numSensors,numElements);
Rz = zeros(numSensors,numElements);

coord = [    0.0000000000000000  0.0000000000000000  0.0000000000000000
    0.0000000000000000  0.0000000000000000  1.0000000000000000
    0.0000000000000000  1.0000000000000000  0.0000000000000000
    1.0000000000000000  0.0000000000000000  0.0000000000000000
    0.0000000000000000  0.0000000000000000  0.7500000000000000
    0.0000000000000000  0.0000000000000000  0.2500000000000000
    0.0000000000000000  0.7500000000000000  0.0000000000000000
    0.0000000000000000  0.7500000000000000  0.2500000000000000
    0.0000000000000000  0.2500000000000000  0.0000000000000000
    0.0000000000000000  0.2500000000000000  0.7500000000000000
    0.7500000000000000  0.0000000000000000  0.0000000000000000
    0.7500000000000000  0.0000000000000000  0.2500000000000000
    0.7500000000000000  0.2500000000000000  0.0000000000000000
    0.2500000000000000  0.0000000000000000  0.0000000000000000
    0.2500000000000000  0.0000000000000000  0.7500000000000000
    0.2500000000000000  0.7500000000000000  0.0000000000000000
    0.5000000000000000  0.5000000000000000  0.0000000000000000
    0.5000000000000000  0.0000000000000000  0.5000000000000000
    0.5000000000000000  0.0000000000000000  0.0000000000000000
    0.0000000000000000  0.5000000000000000  0.5000000000000000
    0.0000000000000000  0.5000000000000000  0.0000000000000000
    0.0000000000000000  0.0000000000000000  0.5000000000000000
    0.2500000000000000  0.2500000000000000  0.0000000000000000
    0.2500000000000000  0.2500000000000000  0.5000000000000000
    0.2500000000000000  0.0000000000000000  0.2500000000000000
    0.2500000000000000  0.0000000000000000  0.5000000000000000
    0.2500000000000000  0.5000000000000000  0.2500000000000000
    0.2500000000000000  0.5000000000000000  0.0000000000000000
    0.0000000000000000  0.2500000000000000  0.2500000000000000
    0.0000000000000000  0.2500000000000000  0.5000000000000000
    0.0000000000000000  0.5000000000000000  0.2500000000000000
    0.5000000000000000  0.2500000000000000  0.2500000000000000
    0.5000000000000000  0.2500000000000000  0.0000000000000000
    0.5000000000000000  0.0000000000000000  0.2500000000000000
    0.2500000000000000  0.2500000000000000  0.2500000000000000];

weights = [  -0.0119047619047619
    -0.0119047619047619
    -0.0119047619047619
    -0.0119047619047619
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    -0.0285714285714286
    -0.0285714285714286
    -0.0285714285714286
    -0.0285714285714286
    -0.0285714285714286
    -0.0285714285714286
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.3047619047619048];

for m = 1:numSensors
    rm = sensor_locations(m,:);
    fun = @(x,y,z) (rm - [x,y,z])./ (sum((rm - [x, y, z]).^2)^1.5);
    for k = 1:numElements

        %Find the vertices of the tetrahedron
        v = fmdl.nodes(fmdl.elems(k,:),:);

        J = [(v(2,:)-v(1,:))',(v(3,:)-v(1,:))',(v(4,:)-v(1,:))'];

        detJ = abs(det(J));

        R = 0;
        for i = 1:length(weights)

            r = coord(i,1);
            s = coord(i,2);
            t = coord(i,3);

            xi = v(1,:)' + J * [r; s; t];

            R = R + weights(i)*fun(xi(1),xi(2),xi(3));
        end

        R =  (detJ / 6) * R;
        Rx(m,k) = R(1);
        Ry(m,k) = R(2);
        Rz(m,k) = R(3);
    end
end

fmdl.R.Rx = Rx;
fmdl.R.Ry = Ry;
fmdl.R.Rz = Rz;

end

function dRdp = dRmkj_dp(fmdl,j,p)
coord = [    0.0000000000000000  0.0000000000000000  0.0000000000000000
    0.0000000000000000  0.0000000000000000  1.0000000000000000
    0.0000000000000000  1.0000000000000000  0.0000000000000000
    1.0000000000000000  0.0000000000000000  0.0000000000000000
    0.0000000000000000  0.0000000000000000  0.7500000000000000
    0.0000000000000000  0.0000000000000000  0.2500000000000000
    0.0000000000000000  0.7500000000000000  0.0000000000000000
    0.0000000000000000  0.7500000000000000  0.2500000000000000
    0.0000000000000000  0.2500000000000000  0.0000000000000000
    0.0000000000000000  0.2500000000000000  0.7500000000000000
    0.7500000000000000  0.0000000000000000  0.0000000000000000
    0.7500000000000000  0.0000000000000000  0.2500000000000000
    0.7500000000000000  0.2500000000000000  0.0000000000000000
    0.2500000000000000  0.0000000000000000  0.0000000000000000
    0.2500000000000000  0.0000000000000000  0.7500000000000000
    0.2500000000000000  0.7500000000000000  0.0000000000000000
    0.5000000000000000  0.5000000000000000  0.0000000000000000
    0.5000000000000000  0.0000000000000000  0.5000000000000000
    0.5000000000000000  0.0000000000000000  0.0000000000000000
    0.0000000000000000  0.5000000000000000  0.5000000000000000
    0.0000000000000000  0.5000000000000000  0.0000000000000000
    0.0000000000000000  0.0000000000000000  0.5000000000000000
    0.2500000000000000  0.2500000000000000  0.0000000000000000
    0.2500000000000000  0.2500000000000000  0.5000000000000000
    0.2500000000000000  0.0000000000000000  0.2500000000000000
    0.2500000000000000  0.0000000000000000  0.5000000000000000
    0.2500000000000000  0.5000000000000000  0.2500000000000000
    0.2500000000000000  0.5000000000000000  0.0000000000000000
    0.0000000000000000  0.2500000000000000  0.2500000000000000
    0.0000000000000000  0.2500000000000000  0.5000000000000000
    0.0000000000000000  0.5000000000000000  0.2500000000000000
    0.5000000000000000  0.2500000000000000  0.2500000000000000
    0.5000000000000000  0.2500000000000000  0.0000000000000000
    0.5000000000000000  0.0000000000000000  0.2500000000000000
    0.2500000000000000  0.2500000000000000  0.2500000000000000];

weights = [  -0.0119047619047619
    -0.0119047619047619
    -0.0119047619047619
    -0.0119047619047619
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    -0.0285714285714286
    -0.0285714285714286
    -0.0285714285714286
    -0.0285714285714286
    -0.0285714285714286
    -0.0285714285714286
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.3047619047619048];

numElements = size(fmdl.elems,1);
numSensors  = numel(fmdl.sensors);
numQuadPts  = length(weights);

% --- Precompute element vertices ---
% V: 4 x 3 x numElements
% Stack vertices: 4 x numElements x 3
V = reshape(fmdl.nodes(fmdl.elems',:), [4, numElements, 3]);

v1 = squeeze(V(1,:,:));  % numElements x 3
v2 = squeeze(V(2,:,:));
v3 = squeeze(V(3,:,:));
v4 = squeeze(V(4,:,:));

% Initialize J
J = zeros(3,3,numElements);  % 3 x 3 x numElements

% Fill J using transpose of row vectors
J(:,1,:) = permute(v2 - v1, [2 3 1]);  % 3 x 1 x numElements
J(:,2,:) = permute(v3 - v1, [2 3 1]);
J(:,3,:) = permute(v4 - v1, [2 3 1]);

% Determinants of all J
detJ = zeros(1,numElements);
for k = 1:numElements
    detJ(k) = abs(det(J(:,:,k)));
end

% --- Compute dRdp for all sensors ---
dRdp = zeros(numSensors,numElements);

for m = 1:numSensors
    rm = fmdl.sensors(m).position;  % 1 x 3
    dm = @(X) rm - X;                % vectorized over X: N x 3

    % Evaluate quadrature points in all elements
    % xi: numQuadPts x 3 x numElements
    xi = reshape(v1', [3,1,numElements]) + ...
        J(:,1,:).*reshape(coord(:,1),[1,numQuadPts,1]) + ...
        J(:,2,:).*reshape(coord(:,2),[1,numQuadPts,1]) + ...
        J(:,3,:).*reshape(coord(:,3),[1,numQuadPts,1]); % 3 x numQuadPts x numElements

    xi = permute(xi,[2 1 3]); % numQuadPts x 3 x numElements

    % dm_vec = rm - xi, same size
    dm_vec = rm - xi; % numQuadPts x 3 x numElements

    % Compute fun values at all quadrature points
    % ((j==p)*norm(dm)^2 - 3*dm_j*dm_p) / norm(dm)^5
    dm_norm2 = sum(dm_vec.^2,2);           % numQuadPts x 1 x numElements
    dm_j    = dm_vec(:,j,:);               % numQuadPts x 1 x numElements
    dm_p    = dm_vec(:,p,:);

    funvals = ((j==p)*dm_norm2 - 3*dm_j.*dm_p) ./ (dm_norm2.^(5/2)); % numQuadPts x 1 x numElements

    % Integrate over quadrature points using weights
    dRdp(m,:) = squeeze(sum(funvals .* reshape(weights, [numQuadPts,1,1]),1))' .* (detJ/6);
end


end

function [dlambda,dR1,dR2] = compute_dlambda_ds(img,dim,p)

num_sensors = numel(img.fwd_model.sensors);
n_nodes = size(img.fwd_model.nodes,1);

dlambda = zeros(n_nodes,num_sensors);

mu_factor = img.fwd_model.mu0/(4*pi);

G = img.fwd_model.G;

sigma = img.elem_data;
A_matrix = M(img,sigma);

% Jacobi preconditioner - matrix free
d = sqrt(diag(A_matrix));        % vector of diagonal entries

Mfun = @(x) x ./ d;              % left preconditioner  M^{-1} x
Nfun = @(x) x ./ d;              % right preconditioner

switch dim
    case 1
        dR1 = dRmkj_dp(img.fwd_model, 3, p);
        G1 = G.Gy;

        dR2 = dRmkj_dp(img.fwd_model, 2, p);
        G2 = G.Gz;
    case 2
        dR1 = dRmkj_dp(img.fwd_model, 1, p);
        dR2 = dRmkj_dp(img.fwd_model, 3, p);

        G1 = G.Gz;
        G2 = G.Gx;
    case 3
        dR1 = dRmkj_dp(img.fwd_model, 2, p);
        dR2 = dRmkj_dp(img.fwd_model, 1, p);

        G1 = G.Gx;
        G2 = G.Gy;
    otherwise
        error('Invalid dimension');
end


% Multiply elementwise by sigma (diagonal) without creating sparse diagonal matrix
% rhs = mu_factor*( dR1 * Sigma * G1 -  dR2* Sigma * G2);
rhs = mu_factor*( (dR1 .* sigma.')*G1 - (dR2 .* sigma.')*G2 );


parfor m = 1:num_sensors
    [dlambda(:,m),~,~] = pcg(A_matrix,rhs(m,:)',1e-9,numel(sigma),Mfun,Nfun);
end

end

function dB = compute_dB_ds(img,dim,p)

num_sensors = numel(img.fwd_model.sensors);
n_nodes = size(img.fwd_model.nodes,1);
n_elem = size(img.fwd_model.elems,1);
num_stim = numel(img.fwd_model.stimulation);

mu_factor = img.fwd_model.mu0/(4*pi);

R = img.fwd_model.R;
G = img.fwd_model.G;

% sigma = img.elem_data;

% Compute EIT forward solution for each current injection pattern
u = fwd_solve(img);
u = u.volt;

switch dim
    case 1
        dRzdp = dRmkj_dp(img.fwd_model, 3, p);

        dRydp = dRmkj_dp(img.fwd_model, 2, p);
    otherwise
        error('Here');
end

Sigma = spdiags(img.elem_data(:), 0, length(img.elem_data), length(img.elem_data));

dB = zeros(num_sensors,num_stim);

for j = 1:num_stim
    % Gx_times_u = G.Gx*u(:,j);

    Gy_times_u = G.Gy*u(:,j);
    Gz_times_u = G.Gz*u(:,j);

    dB(:,j) = -mu_factor*( -dRydp*Sigma*Gz_times_u + dRzdp*Sigma*Gy_times_u);
end

end

function [dJ,dJ1,dJ2] = compute_dJ_ds(img,dim,p)
% Computes derivative of the Jacobian w.r.t. sensor positions
% Vectorized over elements and stimulations, loop only over sensors

num_sensors = numel(img.fwd_model.sensors);
n_elem      = size(img.fwd_model.elems,1);
num_stim    = numel(img.fwd_model.stimulation);

mu_factor = img.fwd_model.mu0/(4*pi);

G = img.fwd_model.G;

% Element volumes
elemV = img.fwd_model.elem_volume(:);  % n_elem x 1

% Compute EIT forward solution for each current injection pattern
u = fwd_solve(img);
u = u.volt;

% G*u for all stim, size: n_elem x num_stim
GxU = G.Gx*u;
GyU = G.Gy*u;
GzU = G.Gz*u;

% Compute derivative of lambda w.r.t sensor positions (output dR1 and dR2
% to avoid recomputing)
[dlambda,dR1dp,dR2dp] = compute_dlambda_ds(img,dim,p);

% Compute derivatives of R w.r.t sensor positions
switch dim
    case 1
        % dR1dp = dRmkj_dp(img.fwd_model, 3, p);
        % dR2dp = dRmkj_dp(img.fwd_model, 2, p);

        G1U = GyU;
        G2U = GzU;
    case 2
        % dR1dp = dRmkj_dp(img.fwd_model, 1, p);
        % dR2dp = dRmkj_dp(img.fwd_model, 3, p);

        G1U = GzU;
        G2U = GxU;
    case 3
        % dR1dp = dRmkj_dp(img.fwd_model, 2, p);
        % dR2dp = dRmkj_dp(img.fwd_model, 1, p);

        G1U = GxU;
        G2U = GyU;
    otherwise
        error('Dimension not supported.');
end


% Preallocate
dJ1 = zeros(num_sensors,num_stim,n_elem);
dJ2 = zeros(num_sensors,num_stim,n_elem);



%% Loop only over sensors
for m = 1:num_sensors
    %% --- dJ1: contribution from dlambda ---
    % dlambda * G matrices, size: n_elem x 1
    dlGx = G.Gx*dlambda(:,m);
    dlGy = G.Gy*dlambda(:,m);
    dlGz = G.Gz*dlambda(:,m);

    % Elementwise multiplication and sum over components
    % tmp: n_elem x num_stim
    tmp = dlGx .* GxU + dlGy .* GyU + dlGz .* GzU;

    % Multiply by element volumes, permute to [num_stim x n_elem]
    dJ1(m,:,:) = tmp.' .* elemV(:).';

    %% --- dJ2: contribution from dR/dp ---
    % dRydp, dRzdp: 1 x n_elem
    % Multiply by G*u per stimulation, permute to [num_stim x n_elem]
    dJ2(m,:,:) = -mu_factor * ( ...
        dR2dp(m,:) .* G2U.' - dR1dp(m,:) .* G1U.' ...
        );
end

%% Total derivative
dJ = dJ1 - dJ2;

end

function J = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,select_sensor_axis)

mu0 = img.fwd_model.mu0;

n_nodes =  size(img.fwd_model.nodes,1);
n_elem = size(img.fwd_model.elems,1);

num_stim = numel(img.fwd_model.stimulation);
num_sensors = numel(img.fwd_model.sensors);

% Compute Gamma matrices
img = compute_gamma_matrices_local(img);

R = img.fwd_model.R;
G = img.fwd_model.G;

switch select_sensor_axis
    case 1
        Gamma = img.Gamma1;
        R1 = R.Rz.';
        R2 = R.Ry.';
    case 2
        Gamma = img.Gamma2;
        R1 = R.Rx.';
        R2 = R.Rz.';
    case 3
        Gamma = img.Gamma3;
        R1 = R.Ry.';
        R2 = R.Rx.';
    otherwise
        error('here')
end

% Compute EIT forward solution for each current injection pattern
u = fwd_solve(img);
u = u.volt;

% Solve the adjoint problem for each sensor to get lambda vectors
lambda = zeros(n_nodes,num_sensors);

A_matrix = A(img.elem_data);

% Jacobi preconditioner - matrix free
d = sqrt(diag(A_matrix));        % vector of diagonal entries

Mfun = @(x) x ./ d;              % left preconditioner  M^{-1} x
Nfun = @(x) x ./ d;              % right preconditioner

GammaT = Gamma.';

% Incomplete Cholesky factorization preconditioner seems to be a bit faster
% than Jacobi preconditioner. However, it breaks down when the
% conductivities become negative
% R = ichol(A_matrix);
% Rt = R';

num_elements = numel(img.elem_data);
for m = 1:num_sensors
    [lambda(:,m),~,~] = pcg(A_matrix,-GammaT(:,m),1e-6,num_elements,Mfun,Nfun);
end

Gx_times_lambda = G.Gx*lambda;
Gy_times_lambda = G.Gy*lambda;
Gz_times_lambda = G.Gz*lambda;

Gx_times_u = G.Gx*u;
Gy_times_u = G.Gy*u;
Gz_times_u = G.Gz*u;

mu_factor = mu0/(4*pi);

elemV = img.fwd_model.elem_volume(:);      % [numElems × 1]

% Expand elem_volume to cover stim × sensor
elemV = reshape(elemV, [n_elem 1 1]);
% Later this will broadcast to [numElems × numStim × numSensors]

% Expand lambda and R terms to 3D
GxL = reshape(Gx_times_lambda, [n_elem 1 num_sensors]); % [: × 1 × numSensors]
GyL = reshape(Gy_times_lambda, [n_elem 1 num_sensors]);
GzL = reshape(Gz_times_lambda, [n_elem 1 num_sensors]);

% Expand u-terms to 3D
GxU = reshape(Gx_times_u, [n_elem num_stim 1]); % [: × numStim × 1]
GyU = reshape(Gy_times_u, [n_elem num_stim 1]);
GzU = reshape(Gz_times_u, [n_elem num_stim 1]);

% Compute all dfdx for all sensors+stim
dfdx = elemV .* ( ...
    GxL.*GxU + ...
    GyL.*GyU + ...
    GzL.*GzU );

% Compute all dfdp (also 3D)
Rx_ = reshape(R.Rx.', [n_elem 1 num_sensors]);
Ry_ = reshape(R.Ry.', [n_elem 1 num_sensors]);
Rz_ = reshape(R.Rz.', [n_elem 1 num_sensors]);

% These are the derivatives with respect to sigma of the C components of the Gamma matrix,
dCxdp = ( -Rz_.*GyU + Ry_.*GzU );
dCydp = ( -Rx_.*GzU + Rz_.*GxU );
dCzdp = ( -Ry_.*GxU + Rx_.*GyU );

% The g matrix does not depend on sigma.
g = zeros(num_sensors,3,3);
for m = 1:numel(img.fwd_model.sensors)
    g(m,:,:) = [...
        img.fwd_model.sensors(m).axes.axis1;
        img.fwd_model.sensors(m).axes.axis2;
        img.fwd_model.sensors(m).axes.axis3];
end

% g: [num_sensors × 3 × 3]
gx = reshape(g(:,select_sensor_axis,1), [1 1 num_sensors]);
gy = reshape(g(:,select_sensor_axis,2), [1 1 num_sensors]);
gz = reshape(g(:,select_sensor_axis,3), [1 1 num_sensors]);

dfdp = mu_factor*(...
    gx.*dCxdp +...
    gy.*dCydp +...
    gz.*dCzdp);

dfd = dfdx + dfdp;   % size: [numElems × numStim × numSensors]

% Now reshape to match J(ids,:)
% permute to [numSensors × numStim × numElems]
dfd = permute(dfd, [3 2 1]);

% collapse first 2 dims → [numSensors*numStim × numElems]
J = reshape(dfd, num_sensors*num_stim, n_elem);

return
end




function cost = compute_cost_function(img,sensor_locations,Gamma_prior,Gamma_noise,A,dim)

n_sensors = numel(img.fwd_model.sensors);
n_stim = numel(img.fwd_model.stimulation);
n_elem = size(img.fwd_model.elems,1);

% Assign sensor locations
img = assign_sensor_locations(img,sensor_locations);

% Compute the jacobian at current sensor locations
J = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,dim);

% Define the inverse posterior covariance matrix
H = J.'*inv(Gamma_noise)*J+inv(Gamma_prior);

cost = trace(inv(H));
end


function dphidp = compute_cost_function_gradient(...
    img,sensor_locations,Gamma_prior,Gamma_noise,A,dim,p)

n_sensors = numel(img.fwd_model.sensors);
n_stim = numel(img.fwd_model.stimulation);
n_elem = size(img.fwd_model.elems,1);

% Assign sensor locations
img = assign_sensor_locations(img,sensor_locations);

% Compute the jacobian at current sensor locations
J = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,dim);

% Compute the Jacobian derivative w.r.t sensor positions
dJds = compute_dJ_ds(img,dim,p);

% Define the inverse posterior covariance matrix
H = J.'*inv(Gamma_noise)*J+inv(Gamma_prior);

% Intermediate variable W
W =  inv(Gamma_noise)*J*(inv(H)^2); %[n_stim*n_sensors,n_elem]

% Compute the derivative of the cost w.r.t. the p-coordinate of
% the sensor
dphidp = zeros(1,n_sensors);
for m = 1:n_sensors
    temp = 0.0;
    for l = 1:n_stim
        id = m+(l-1)*n_sensors;
        temp = temp -2*sum(W(id,:) .* squeeze(dJds(m,l,:)).'); % Frobenius inner product
    end
    dphidp(m) = temp;
end

end
%% Helper functions

%function: Assign sensor locations
function img = assign_sensor_locations(img,sensor_locations)
assert(numel(img.fwd_model.sensors) == size(sensor_locations,1));
for id = 1: numel(img.fwd_model.sensors)
    img.fwd_model.sensors(id).position = sensor_locations(id,:);
end
end

function sensor_locations = fetch_sensor_locations(img)
n_sensors = numel(img.fwd_model.sensors);
sensor_locations = zeros(n_sensors,3);

for m = 1: numel(img.fwd_model.sensors)
    sensor_locations(m,:) = img.fwd_model.sensors(m).position;
end

end

function x = sensor_locations_to_vector(sensor_locations)

n_sensors = size(sensor_locations,1);

x = zeros(n_sensors*3,1);

for m = 1:n_sensors
    for p = 1:3
        id = m + (p-1)*n_sensors;
        x(id) = sensor_locations(m,p);
    end
end

end

function sensor_locations = vector_to_sensor_locations(x)

assert(mod(numel(x),3) == 0);

n_sensors = numel(x)/3;

sensor_locations = zeros(n_sensors,3);
for m = 1:n_sensors
    for p = 1:3
        id = m + (p-1)*n_sensors;
        sensor_locations(m,p) = x(id);
    end
end
end

function dJ = compute_dJ_ds_complex(img,A,sensor_locations_0,delta,dim,p)
n_stim = numel(img.fwd_model.stimulation);
n_sensors = numel(img.fwd_model.sensors);
n_elem = size(img.fwd_model.elems,1);

dJ = zeros(n_sensors,n_stim,n_elem);

for m = 1:numel(img.fwd_model.sensors)
    fprintf('Done iteration %i of %i\n',m,numel(img.fwd_model.sensors))

    sensor_locations_perturbed = sensor_locations_0;
    sensor_locations_perturbed(m,p) = sensor_locations_0(m,p) + 1i*delta;
    img_perturbed = assign_sensor_locations(img,sensor_locations_perturbed);

    jacobian_perturbed = calc_jacobian_1axis_direct_fully_vectorized_local(img_perturbed,A,dim);

    for l = 1:n_stim
        id = m+n_sensors*(l-1);
        dJ(m,l,:) = imag(jacobian_perturbed(id,:)) / delta;
    end
end
end

function dJ = compute_dJ_ds_central_differences(img,A,sensor_locations_0,delta,dim,p)

n_stim = numel(img.fwd_model.stimulation);
n_sensors = numel(img.fwd_model.sensors);
n_elem = size(img.fwd_model.elems,1);

dJ = zeros(n_sensors,n_stim,n_elem);

img_base = img;   % freeze once outside


for m = 1:numel(img.fwd_model.sensors)
    fprintf('Done iteration %i of %i\n',m,numel(img.fwd_model.sensors))

    img_plus  = img_base;
    img_minus = img_base;

    sensor_locations_perturbed = sensor_locations_0;
    sensor_locations_perturbed(m,p) = sensor_locations_0(m,p) + delta;
    img_plus_perturbed = assign_sensor_locations(img_plus,sensor_locations_perturbed);

    jacobian_x_perturbed_p_plus = calc_jacobian_1axis_direct_fully_vectorized_local(img_plus_perturbed,A,dim);

    sensor_locations_perturbed = sensor_locations_0;
    sensor_locations_perturbed(m,p) = sensor_locations_0(m,p) - delta;
    img_minus_perturbed = assign_sensor_locations(img_minus,sensor_locations_perturbed);

    jacobian_x_perturbed_p_minus = calc_jacobian_1axis_direct_fully_vectorized_local(img_minus_perturbed,A,dim);

    for l = 1:n_stim
        id = m+n_sensors*(l-1);
        dJ(m,l,:) = (jacobian_x_perturbed_p_plus(id,:)-jacobian_x_perturbed_p_minus(id,:))/(2*delta);
    end
end
end