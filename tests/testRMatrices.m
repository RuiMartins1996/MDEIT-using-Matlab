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


            Gamma_prior_inv = inv(Gamma_prior);
            Gamma_noise_inv = inv(Gamma_noise);

            %% Compare analytical and numerical computation of gradient
            A = @(x) M(img,x);
            delta = 1e-5;

            errorThresh = testCase.testParameters.errorThresh;

            function dphidp_numerical = compute_cost_function_gradient_central_differences(img,A,dim,p,Gamma_prior_inv,Gamma_noise_inv,delta)

                n_sensors = numel(img.fwd_model.sensors);
                dphidp_numerical = zeros(1,n_sensors);

                sensor_locations_0_local = fetch_sensor_locations(img);

                for sensor_id = 1:n_sensors
                    sensor_locations_local = sensor_locations_0_local;
                    sensor_locations_local(sensor_id,p) = sensor_locations_local(sensor_id,p) + delta;

                    % assign sensor locations is done inside
                    % compute_cost_function
                    cost_m_plus = compute_cost_function_a_opt(img,sensor_locations_local,Gamma_prior_inv,Gamma_noise_inv,A,dim);

                    sensor_locations_local = sensor_locations_0_local;
                    sensor_locations_local(sensor_id,p) = sensor_locations_local(sensor_id,p) - delta;

                    cost_m_minus = compute_cost_function_a_opt(img,sensor_locations_local,Gamma_prior_inv,Gamma_noise_inv,A,dim);

                    dphidp_numerical(sensor_id) = (cost_m_plus-cost_m_minus)/(2*delta);
                end
            end

            for dim = 1:3
                for p = 1:3

                    %Compute the gradient of the cost function analytically
                    dphidp = compute_cost_function_gradient_a_opt(...
                        img,sensorLocations,Gamma_prior_inv,Gamma_noise_inv ,A,dim,p);

                    %Compute the gradient of the cost function numerically with central-differences
                    dphidp_numerical = compute_cost_function_gradient_central_differences(img,A,dim,p,Gamma_prior_inv,Gamma_noise_inv,delta);

                    err = norm(dphidp_numerical(:) - dphidp(:)) / norm(dphidp(:));
                    testCase.verifyTrue(all(err<errorThresh),sprintf('Analytical and numerical difference is more than %d (%%)\n',errorThresh));
                end
            end

        end

        function testCostFunctionGradientWithRespectToCylindricalCoordinates(testCase)
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

            Gamma_prior_inv = inv(Gamma_prior);
            Gamma_noise_inv = inv(Gamma_noise);

            %% Compare analytical and numerical computation of gradient
            A = @(x) M(img,x);
            delta = 1e-9;

            errorThresh = testCase.testParameters.errorThresh;

            function [dphidr_numerical,dphidtheta_numerical,dphidz_numerical] = ...
                    compute_cost_function_gradient_central_differences(img,A,dim,Gamma_prior_inv,Gamma_noise_inv,delta)

                n_sensors = numel(img.fwd_model.sensors);
                dphidr_numerical = zeros(1,n_sensors);
                dphidtheta_numerical = zeros(1,n_sensors);
                dphidz_numerical = zeros(1,n_sensors);

                sensor_locations_0_local = fetch_sensor_locations(img);

                % Derivative w.r.t. r
                for sensor_id = 1:n_sensors
                    sensor_locations_local = sensor_locations_0_local;

                    rm =  sqrt(sensor_locations_local(sensor_id,1)^2+sensor_locations_local(sensor_id,2)^2);
                    thetam = atan2(sensor_locations_local(sensor_id,2),sensor_locations_local(sensor_id,1));
                    zm =  sensor_locations_local(sensor_id,3);

                    % Positive delta
                    sensor_locations_local(sensor_id,:) = ...
                        [(rm+delta)*cos(thetam),(rm+delta)*sin(thetam),zm];

                    cost_m_plus = compute_cost_function_a_opt(img,sensor_locations_local,Gamma_prior_inv,Gamma_noise_inv,A,dim);

                    % Negative delta
                    sensor_locations_local(sensor_id,:) = ...
                        [(rm-delta)*cos(thetam),(rm-delta)*sin(thetam),zm];
                    cost_m_minus = compute_cost_function_a_opt(img,sensor_locations_local,Gamma_prior_inv,Gamma_noise_inv,A,dim);

                    dphidr_numerical(sensor_id) = (cost_m_plus-cost_m_minus)/(2*delta);
                end

                % Derivative w.r.t. theta
                for sensor_id = 1:n_sensors
                    sensor_locations_local = sensor_locations_0_local;

                    rm =  sqrt(sensor_locations_local(sensor_id,1)^2+sensor_locations_local(sensor_id,2)^2);
                    thetam = atan2(sensor_locations_local(sensor_id,2),sensor_locations_local(sensor_id,1));
                    zm =  sensor_locations_local(sensor_id,3);

                    % Positive delta
                    sensor_locations_local(sensor_id,:) = ...
                        [rm*cos(thetam+delta),rm*sin(thetam+delta),zm];

                    cost_m_plus = compute_cost_function_a_opt(img,sensor_locations_local,Gamma_prior_inv,Gamma_noise_inv,A,dim);

                    % Negative delta
                    sensor_locations_local(sensor_id,:) = ...
                        [rm*cos(thetam-delta),rm*sin(thetam-delta),zm];
                    cost_m_minus = compute_cost_function_a_opt(img,sensor_locations_local,Gamma_prior_inv,Gamma_noise_inv,A,dim);

                    dphidtheta_numerical(sensor_id) = (cost_m_plus-cost_m_minus)/(2*delta);
                end

                % Derivative w.r.t. z
                for sensor_id = 1:n_sensors
                    sensor_locations_local = sensor_locations_0_local;

                    rm =  sqrt(sensor_locations_local(sensor_id,1)^2+sensor_locations_local(sensor_id,2)^2);
                    thetam = atan2(sensor_locations_local(sensor_id,2),sensor_locations_local(sensor_id,1));
                    zm =  sensor_locations_local(sensor_id,3);

                    % Positive delta
                    sensor_locations_local(sensor_id,:) = ...
                        [rm*cos(thetam),rm*sin(thetam),zm+delta];

                    cost_m_plus = compute_cost_function_a_opt(img,sensor_locations_local,Gamma_prior_inv,Gamma_noise_inv,A,dim);

                    % Negative delta
                    sensor_locations_local(sensor_id,:) = ...
                        [rm*cos(thetam),rm*sin(thetam),zm-delta];
                    cost_m_minus = compute_cost_function_a_opt(img,sensor_locations_local,Gamma_prior_inv,Gamma_noise_inv,A,dim);

                    dphidz_numerical(sensor_id) = (cost_m_plus-cost_m_minus)/(2*delta);
                end
            end

            function [dphidr_complex,dphidtheta_complex,dphidz_complex] = ...
                    compust_cost_function_gradient_complex(img,A,dim,Gamma_prior_inv,Gamma_noise_inv,delta)

                n_stim = numel(img.fwd_model.stimulation);
                n_sensors = numel(img.fwd_model.sensors);
                n_elem = size(img.fwd_model.elems,1);

                dphidr_complex = zeros(1,n_sensors);
                dphidtheta_complex = zeros(1,n_sensors);
                dphidz_complex = zeros(1,n_sensors);

                sensor_locations_0_local = fetch_sensor_locations(img);

                % Derivative with respect to r
                for m = 1:numel(img.fwd_model.sensors)
                    sensor_locations_local = sensor_locations_0_local;

                    rm =  sqrt(sensor_locations_local(m,1)^2+sensor_locations_local(m,2)^2);
                    thetam = atan2(sensor_locations_local(m,2),sensor_locations_local(m,1));
                    zm =  sensor_locations_local(m,3);

                    sensor_locations_local(m,:) = ...
                        [(rm+1i*delta)*cos(thetam),(rm+1i*delta)*sin(thetam),zm];

                    cost_perturbed = compute_cost_function_a_opt(img,sensor_locations_local,Gamma_prior_inv,Gamma_noise_inv,A,dim);

                    dphidr_complex(m) = imag(cost_perturbed) / delta;
                end

                % Derivative with respect to theta
                for m = 1:numel(img.fwd_model.sensors)
                    sensor_locations_local = sensor_locations_0_local;

                    rm =  sqrt(sensor_locations_local(m,1)^2+sensor_locations_local(m,2)^2);
                    thetam = atan2(sensor_locations_local(m,2),sensor_locations_local(m,1));
                    zm =  sensor_locations_local(m,3);

                    sensor_locations_local(m,:) = ...
                        [rm*cos(thetam+1i*delta),rm*sin(thetam+1i*delta),zm];

                    cost_perturbed = compute_cost_function_a_opt(img,sensor_locations_local,Gamma_prior_inv,Gamma_noise_inv,A,dim);

                    dphidtheta_complex(m) = imag(cost_perturbed) / delta;
                end

                % Derivative with respect toz
                for m = 1:numel(img.fwd_model.sensors)
                    sensor_locations_local = sensor_locations_0_local;

                    rm =  sqrt(sensor_locations_local(m,1)^2+sensor_locations_local(m,2)^2);
                    thetam = atan2(sensor_locations_local(m,2),sensor_locations_local(m,1));
                    zm =  sensor_locations_local(m,3);

                    sensor_locations_local(m,:) = ...
                        [rm*cos(thetam),rm*sin(thetam),zm+1i*delta];

                    cost_perturbed = compute_cost_function_a_opt(img,sensor_locations_local,Gamma_prior_inv,Gamma_noise_inv,A,dim);

                    dphidz_complex(m) = imag(cost_perturbed) / delta;
                end
            end

            function [dphidr,dphidtheta,dphidz] = compute_cost_function_gradient(img,A,dim,Gamma_prior_inv,Gamma_noise_inv)

                n_sensors = numel(img.fwd_model.sensors);

                sensor_locations = fetch_sensor_locations(img);

                rm = zeros(n_sensors,1);
                thetam = zeros(n_sensors,1);
                zm = zeros(n_sensors,1);

                for m = 1:n_sensors
                    rm(m) = sqrt(sensor_locations(m,1)^2+sensor_locations(m,2)^2);
                    thetam(m) = atan2(sensor_locations(m,2),sensor_locations(m,1));
                    zm(m) = sensor_locations(m,3);
                end


                dphidp = compute_cost_function_gradient_a_opt_optimized(...
                    img,sensor_locations,Gamma_prior_inv,Gamma_noise_inv,A,dim);

                dphidx = dphidp(1,:);
                dphidy = dphidp(2,:);
                dphidz = dphidp(3,:);

                dphidr = cos(thetam)'.*dphidx+sin(thetam)'.*dphidy;
                dphidtheta = -rm'.*sin(thetam)'.*dphidx + rm'.*cos(thetam)'.*dphidy;

            end

            for dim = 1:3
                %Compute the gradient of the cost function analytically
                [dphidr,dphidtheta,dphidz] = ...
                    compute_cost_function_gradient(img,A,dim,Gamma_prior_inv,Gamma_noise_inv);

                %Compute the gradient of the cost function numerically with central-differences
                [dphidr_numerical,dphidtheta_numerical,dphidz_numerical] = ...
                    compute_cost_function_gradient_central_differences(img,A,dim,Gamma_prior_inv,Gamma_noise_inv,delta);

                %Compute the gradient of the cost function numerically with
                %complex step method
                [dphidr_complex,dphidtheta_complex,dphidz_complex] = ...
                    compust_cost_function_gradient_complex(img,A,dim,Gamma_prior_inv,Gamma_noise_inv,delta);

                err_r_1 = 100*max(abs(dphidr_numerical(:) - dphidr(:))./ abs(dphidr(:)));
                err_r_2 = 100*max(abs(dphidr_complex(:) - dphidr(:))./ abs(dphidr(:)));

                err_theta_1 = 100*max(abs(dphidtheta_numerical(:) - dphidtheta(:))./ abs(dphidtheta(:)));
                err_theta_2 = 100*max(abs(dphidtheta_complex(:) - dphidtheta(:))./ abs(dphidtheta(:)));

                err_z_1 = 100*max(abs(dphidz_numerical(:) - dphidz(:))./ abs(dphidz(:)));
                err_z_2 = 100*max(abs(dphidz_complex(:) - dphidz(:))./ abs(dphidz(:)));

                testCase.verifyTrue(err_r_1<errorThresh | err_r_2<errorThresh,sprintf('Analytical and numerical difference is more than %d (%%)\n',errorThresh));
                testCase.verifyTrue(err_theta_1<errorThresh | err_theta_2<errorThresh,sprintf('Analytical and numerical difference is more than %d (%%)\n',errorThresh));
                testCase.verifyTrue(err_z_1<errorThresh | err_z_2<errorThresh,sprintf('Analytical and numerical difference is more than %d (%%)\n',errorThresh));

            end

        end

        function testOptimizationDifferentInitialConditions(testCase)
            %% Prepare workspace
            % Get the full path of the current script
            fullpath = mfilename('fullpath');
            % Extract just the folder
            script_folder = fileparts(fullpath);
            cd(script_folder);

            % Have to add the functions path manually so prepare_workspace runs
            grandparent_folder =fileparts(fileparts(script_folder));
            addpath(genpath(fullfile(grandparent_folder,'functions')));

            model_folder = prepare_workspace(script_folder);

            rng(1);

            %% Model parameters ( only two sensors, to reduce dimensionality of the problem so we can plot it)
            z0 = 0.0058; %(Ohm m^2) is the contact impedance from the CEM article 58 Ohm cm^2
            l0 = 40e-3; %(m) the tank radius
            I0 = 2.4e-3;%(A) the magnitude of the injected current

            % The derived characteristic units
            V0 = z0*I0/(l0^2); %(V)
            sigma0 = l0/z0; %(S/m)
            J0 = I0/(l0^2);

            model_parameters = create_default_3d_model_parameters(l0, z0, sigma0, I0);

            model_parameters.maxsz = max(model_parameters.height,model_parameters.radius)/8;
            model_parameters.numOfElectrodesPerRing = 4;
            model_parameters.numOfRings = 2;
            model_parameters.numOfSensors = 2;
            model_parameters.sensorRadius = model_parameters.radius*1.5;
            model_parameters.material.type = 'spherical';
            model_parameters.material.position(1) = 0.95*(model_parameters.radius-model_parameters.material.radius);
            model_parameters.material.position(3) = 0.5*model_parameters.height;

            background_conductivity = 3.28e-1/sigma0;
            anomaly_conductivity = background_conductivity/10;

            %% Simulation parameters
            dim = 3;            %if doing 1-axis, which dimmension?

            %Optimizer parameters
            max_iteratons = 50;
            hessian_approximation = 'lbfgs';
            use_parallel = true;

            rmax = 2*model_parameters.radius;
            rmin = model_parameters.radius*1.1;

            zmax = model_parameters.height;
            zmin = 0;

            r0 = 1.5*model_parameters.radius; %radius of cylinder shell to perform optimizaiton on

            alpha = 0; %controls the force of the repulsion!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


            %% Make forward model

            [~,fmdls] = mk_mdeit_model(model_parameters,model_folder);

            fmdl = fmdls{1};

            sensorLocations = zeros(numel(fmdl.sensors),3);
            for i = 1: numel(fmdl.sensors)
                sensorLocations(i,:) = fmdl.sensors(i).position;
            end

            % Make homogeneous image
            imgh = mk_image_mdeit(fmdl,background_conductivity);

            % Add plastic cylinder
            imgi = add_material_properties(imgh, [background_conductivity,anomaly_conductivity]);

            % show_fem(imgi);
            % plot_sensors(imgi);

            n_sensors = numel(imgi.fwd_model.sensors);
            n_stim = numel(imgi.fwd_model.stimulation);
            n_elem = size(imgi.fwd_model.elems,1);
            n_nodes = size(fmdl.nodes,1);

            A = @(x) M(imgi,x);

            %% Initialize

            n_trial = 2;

            R0 = model_parameters.radius*1.5;
            r_new = R0*ones(1,n_sensors);
            theta_new = [0,pi];
            
            
            all_sensor_locations_0 = zeros(n_sensors,3,n_trial);

            % Set multiple initial conditions
            for i = 1:n_trial
                z_new = model_parameters.height*rand(1,n_sensors);
  
                sensor_locations_0 = [(r_new.*cos(theta_new))',(r_new.*sin(theta_new))',z_new'];
                
                all_sensor_locations_0(:,:,i) = sensor_locations_0;
            end

            
            %% Define prior and noise covariance matrices

            B = fwd_solve_mdeit(imgi);
            max_B = max([abs(B.Bx(:));abs(B.By(:));abs(B.Bz(:))]);

            % Set the noise variance with respect to the data magnitude
            noise_std_deviation = max_B/10;
            variance_noise = noise_std_deviation^2;

            Jdim = calc_jacobian_1axis_direct_fully_vectorized_local(imgi,A,dim);

            coeff = 50;
            variance_prior = coeff*variance_noise/eigs(Jdim'*Jdim,1);

            % Check how many eigenvectors of J'J are in the data-dominated
            % regime, as opposed to the prior-dominated regime
            d = sum(eigs(Jdim'*Jdim,n_elem).*variance_prior/variance_noise>1);

            fprintf('# lambda_i*alpha:beta>1 = %i (%i)\n',d,n_elem);

            % White noise with zero mean and \mu variance
            Gamma_noise = variance_noise.*speye(n_stim*n_sensors,n_stim*n_sensors);
            inv_Gamma_noise = inv(Gamma_noise);

            Gamma_prior = variance_prior.*speye(n_elem,n_elem);
            inv_Gamma_prior = inv(Gamma_prior);

            Jx = calc_jacobian_1axis_direct_fully_vectorized_local(imgi,A,1);
            Jy = calc_jacobian_1axis_direct_fully_vectorized_local(imgi,A,2);
            Jz = calc_jacobian_1axis_direct_fully_vectorized_local(imgi,A,3);

            J_3axis = [Jx;Jy;Jz];

            coeff = 50;
            variance_prior_3_axis = coeff*variance_noise/eigs(J_3axis'*J_3axis,1);

            d_3axis = sum(eigs(J_3axis'*J_3axis,n_elem).*variance_prior_3_axis/variance_noise>1);

            fprintf('# lambda_i*alpha:beta>1 = %i (%i)\n',d_3axis,n_elem);

            Gamma_prior_3_axis = variance_prior_3_axis.*speye(n_elem,n_elem);
            inv_Gamma_prior_3_axis = inv(Gamma_prior_3_axis);

            Gamma_noise_3_axis = variance_noise.*speye(3*n_stim*n_sensors,3*n_stim*n_sensors);
            inv_Gamma_noise_3_axis = inv(Gamma_noise_3_axis);


            %% Define functions
            jacobian_coordinate_transformation_cylindrical = compute_jacobian_coordinate_transformation_cylindrical(sensor_locations_0);

            function jacobian_coordinate_transformation = compute_jacobian_coordinate_transformation_cylindrical(sensor_locations)

                n_sensors = size(sensor_locations,1);

                jacobian_coordinate_transformation = zeros(3,3,n_sensors);

                rm = sqrt(sensor_locations(:,1).^2+sensor_locations(:,2).^2);
                thetam = atan2(sensor_locations(:,2),sensor_locations(:,1));
                zm = sensor_locations(:,3);

                jacobian_coordinate_transformation_r = zeros(3,n_sensors);
                jacobian_coordinate_transformation_r(1,:) = cos(thetam);
                jacobian_coordinate_transformation_r(2,:) = sin(thetam);
                jacobian_coordinate_transformation_r(3,:) = zeros(1,n_sensors);

                jacobian_coordinate_transformation_theta = zeros(3,n_sensors);
                jacobian_coordinate_transformation_theta(1,:) = -rm'.*sin(thetam)';
                jacobian_coordinate_transformation_theta(2,:) = rm'.*cos(thetam)';
                jacobian_coordinate_transformation_theta(3,:) = zeros(1,n_sensors);

                jacobian_coordinate_transformation_z = zeros(3,n_sensors);
                jacobian_coordinate_transformation_z(1,:) = zeros(1,n_sensors);
                jacobian_coordinate_transformation_z(2,:) = zeros(1,n_sensors);
                jacobian_coordinate_transformation_z(3,:) = ones(1,n_sensors);

                jacobian_coordinate_transformation(1,:,:) = jacobian_coordinate_transformation_r;
                jacobian_coordinate_transformation(2,:,:) = jacobian_coordinate_transformation_theta;
                jacobian_coordinate_transformation(3,:,:) = jacobian_coordinate_transformation_z;
            end

            % Map sensor locations to vector in cartesian coordinates
            function x = sensor_locations_to_vector_cartesian(sensor_locations)

                n_sensors = size(sensor_locations,1);

                x = zeros(3*n_sensors,1);

                x(1:n_sensors) =  sensor_locations(:,1);
                x(n_sensors+1:2*n_sensors) = sensor_locations(:,2);
                x(2*n_sensors+1:3*n_sensors) = sensor_locations(:,3);

            end

            function sensor_locations = vector_to_sensor_locations_cartesian(x)

                assert(mod(numel(x),3)==0);
                n_sensors = numel(x)/3;

                sensor_locations = zeros(n_sensors,3);

                sensor_locations(:,1) = x(1:n_sensors);
                sensor_locations(:,2) = x(n_sensors+1:2*n_sensors);
                sensor_locations(:,3) = x(2*n_sensors+1:3*n_sensors);

            end

            % Map cartesian coordinates to cylindrical coordinates
            function q = x_to_q_cylindrical(x)

                assert(mod(numel(x),3)==0);
                n_sensors = numel(x)/3;

                rm = sqrt(x(1:n_sensors).^2+x(n_sensors+1:2*n_sensors).^2);
                thetam = atan2(x(n_sensors+1:2*n_sensors),x(1:n_sensors));
                zm = x(2*n_sensors+1:3*n_sensors);

                q = [rm(:);thetam(:);zm(:)];
            end

            function x = q_to_x_cylindrical(q)

                assert(mod(numel(q),3)==0);
                n_sensors = numel(q)/3;

                rm = q(1:n_sensors);
                thetam = q(n_sensors+1:2*n_sensors);
                zm = q(2*n_sensors+1:3*n_sensors);

                x = [rm(:).*cos(thetam(:));rm(:).*sin(thetam(:));zm(:)];
            end

            % Map cartesian coordinates to parametrized region (region contained
            % between two cylinders of radius rmin and rmax and height [0,3])
            function q = map_x_to_q_cyl_region(x,rmin,rmax,zmin,zmax)
                assert(mod(numel(x),3)==0);
                n_sensors = numel(x)/3;

                rm = sqrt(x(1:n_sensors).^2+x(n_sensors+1:2*n_sensors).^2);
                thetam = atan2(x(n_sensors+1:2*n_sensors),x(1:n_sensors));
                zm = x(2*n_sensors+1:3*n_sensors);


                xim = log((rm-rmin)./(rmax-rm));
                etam = log((zm-zmin)./(zmax-zm));

                q = [xim(:);thetam(:);etam(:)];
            end

            function x = map_q_to_x_cyl_region(q,rmin,rmax,zmin,zmax)

                sigmoid = @(y) 1./(1+exp(-y));

                assert(mod(numel(q),3)==0);
                n_sensors = numel(q)/3;

                q = q(:);
                xim = q(1:n_sensors);
                thetam = q(n_sensors+1:2*n_sensors);
                etam = q(2*n_sensors+1:3*n_sensors);

                rm = rmin + (rmax-rmin).*sigmoid(xim);
                zm = zmin + (zmax-zmin).*sigmoid(etam);

                x = [rm(:).*cos(thetam(:));rm(:).*sin(thetam(:));zm(:)];
            end

            x_to_q_cyl_region = @(x) map_x_to_q_cyl_region(x,rmin,rmax,zmin,zmax);

            q_to_x_cyl_region = @(q) map_q_to_x_cyl_region(q,rmin,rmax,zmin,zmax);

            % Map cartesian coordinates to parametrized region (cylinder shell at radius r0)
            function q = map_x_to_q_cyl_shell(x,r0,zmin,zmax)
                assert(mod(numel(x),2)==0);
                n_sensors = numel(x)/2;

                thetam = atan2(x(n_sensors+1:2*n_sensors),x(1:n_sensors));
                zm = x(2*n_sensors+1:3*n_sensors);

                etam = log((zm-zmin)./(zmax-zm));

                q = [thetam(:);etam(:)];
            end

            function x = map_q_to_x_cyl_shell(q,r0,zmin,zmax)
                sigmoid = @(y) 1./(1+exp(-y));

                assert(mod(numel(q),2)==0);
                n_sensors = numel(q)/2;

                q = q(:);
                thetam = q(1:n_sensors);
                etam = q(n_sensors+1:2*n_sensors);

                zm = zmin + (zmax-zmin).*sigmoid(etam);

                x = [r0.*cos(thetam(:));r0.*sin(thetam(:));zm(:)];
            end

            x_to_q_cyl_shell = @(x) map_x_to_q_cyl_shell(x,r0,zmin,zmax);

            q_to_x_cyl_shell = @(q) map_q_to_x_cyl_shell(q,r0,zmin,zmax);

            x0 = sensor_locations_to_vector_cartesian(sensor_locations_0);

            % Full map from q coordinates to sensor locations and back for constrained
            % optimization
            vector_to_sensor_locations_con = @(q) vector_to_sensor_locations_cartesian(q_to_x_cylindrical(q));
            sensor_locations_to_vector_con = @(sensor_locations) x_to_q_cylindrical(sensor_locations_to_vector_cartesian(sensor_locations));

            % Full map from q coordinates to sensor locations and back for
            % unconstrained optimization
            vector_to_sensor_locations_unc = @(q) vector_to_sensor_locations_cartesian(q_to_x_cyl_region(q));
            sensor_locations_to_vector_unc = @(sensor_locations) x_to_q_cyl_region(sensor_locations_to_vector_cartesian(sensor_locations));

            % Sanity check
            q = x_to_q_cylindrical(x0);
            x = q_to_x_cylindrical(q);
            assert(norm(x-x0)<1e-5,'Unexpected');

            % Sanity check
            q0 = sensor_locations_to_vector_con(sensor_locations_0);
            sensor_locations_new = vector_to_sensor_locations_con(q0);
            assert(norm(sensor_locations_new-sensor_locations_0)<1e-5,'Unexpected');

            % Sanity check
            sensor_locations = vector_to_sensor_locations_con(q0);
            q_new = sensor_locations_to_vector_con(sensor_locations);
            assert(norm(q0-q_new)<1e-5,'Unexpected');

            % Transform gradient of cost function w.r.t. cartesian coordinates to
            % w.r.t. q coordinates
            function dphidq = dphidp_to_dphidq(sensor_locations,dphidp,jacobian_coordinate_transformation)

                n_sensors = size(sensor_locations,1);

                assert(all(size(jacobian_coordinate_transformation) == [3,3,n_sensors]));

                dphidq = zeros(3,n_sensors);

                for q = 1:3
                    temp = zeros(1,n_sensors);
                    for dim = 1:3
                        temp = temp + squeeze(jacobian_coordinate_transformation(q,dim,:)).'.*dphidp(dim,:);
                    end
                    dphidq(q,:) = temp;
                end

            end


            %% Functions to compute cost and cost gradient
            function out = f(img,q,inv_Gamma_prior,inv_Gamma_noise,A,vector_to_sensor_locations,opt_mode,mode,dim)

                allowed_opt_modes = {'a-opt','d-opt'};
                assert(ismember(opt_mode,allowed_opt_modes));

                allowed_modes = {'mdeit3','mdeit1'};
                assert(ismember(mode,allowed_modes));

                if strcmp(mode,'mdeit3') && nargin<9
                    dim = 'default';
                end

                sensor_locations = vector_to_sensor_locations(q);

                switch mode
                    case 'mdeit1'
                        switch opt_mode
                            case 'a-opt'
                                phi_oed = compute_cost_function_a_opt(img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A,dim);
                            case 'd-opt'
                                phi_oed = compute_cost_function_d_opt(img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A,dim);
                        end
                    case 'mdeit3'
                        switch opt_mode
                            case 'a-opt'
                                phi_oed = compute_cost_function_a_opt_3_axis(img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A);
                            case 'd-opt'
                                phi_oed = compute_cost_function_d_opt_3_axis(img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A);
                        end
                end

                out = phi_oed;
            end

            function dphi = g(img,q,inv_Gamma_prior,inv_Gamma_noise,A,...
                    vector_to_sensor_locations,jacobian_coordinate_transformation,opt_mode,mode,dim)

                n_sensors = numel(img.fwd_model.sensors);

                allowed_opt_modes = {'a-opt','d-opt'};
                assert(ismember(opt_mode,allowed_opt_modes));

                allowed_modes = {'mdeit3','mdeit1'};
                assert(ismember(mode,allowed_modes));

                if strcmp(mode,'mdeit3') && nargin<10
                    dim = 'default';
                end

                sensor_locations = vector_to_sensor_locations(q);

                switch mode
                    case 'mdeit1'
                        switch opt_mode
                            case 'a-opt'
                                dphidp = compute_cost_function_gradient_a_opt_optimized(...
                                    img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A,dim);
                            case 'd-opt'
                                dphidp = compute_cost_function_gradient_d_opt_optimized(...
                                    img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A,dim);
                        end
                    case 'mdeit3'
                        switch opt_mode
                            case 'a-opt'
                                dphidp = compute_cost_function_gradient_a_opt_optimized_3_axis(...
                                    img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A);
                            case 'd-opt'
                                dphidp = compute_cost_function_gradient_d_opt_optimized_3_axis(...
                                    img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A);
                        end
                end

                % Convert cartesian derivatives to other coordinate derivatives (
                % generalized)
                dphidq = dphidp_to_dphidq(sensor_locations,dphidp,jacobian_coordinate_transformation);

                dphi = [dphidq(1,:),dphidq(2,:),dphidq(3,:)];


            end

            f_a_opt_mdeit_dim_con = @(q) f(imgi,q,inv_Gamma_prior,inv_Gamma_noise,A,vector_to_sensor_locations_con,'a-opt','mdeit1',dim);
            g_a_opt_mdeit_dim_con  = @(q) g(imgi,q,inv_Gamma_prior,inv_Gamma_noise,A,...
                vector_to_sensor_locations_con,jacobian_coordinate_transformation_cylindrical,'a-opt','mdeit1',dim);

            f_d_opt_mdeit_dim_con = @(q) f(imgi,q,inv_Gamma_prior,inv_Gamma_noise,A,vector_to_sensor_locations_con,'d-opt','mdeit1',dim);
            g_d_opt_mdeit_dim_con  = @(q) g(imgi,q,inv_Gamma_prior,inv_Gamma_noise,A,...
                vector_to_sensor_locations_con,jacobian_coordinate_transformation_cylindrical,'d-opt','mdeit1',dim);

            f_a_opt_mdeit3_con = @(q) f(imgi,q,inv_Gamma_prior_3_axis,inv_Gamma_noise_3_axis,A,vector_to_sensor_locations_con,'a-opt','mdeit3',dim);
            g_a_opt_mdeit3_con  = @(q) g(imgi,q,inv_Gamma_prior_3_axis,inv_Gamma_noise_3_axis,A,...
                vector_to_sensor_locations_con,jacobian_coordinate_transformation_cylindrical,'a-opt','mdeit3',dim);

            f_d_opt_mdeit3_con = @(q) f(imgi,q,inv_Gamma_prior_3_axis,inv_Gamma_noise_3_axis,A,vector_to_sensor_locations_con,'d-opt','mdeit3',dim);
            g_d_opt_mdeit3_con  = @(q) g(imgi,q,inv_Gamma_prior_3_axis,inv_Gamma_noise_3_axis,A,...
                vector_to_sensor_locations_con,jacobian_coordinate_transformation_cylindrical,'d-opt','mdeit3',dim);

            %% PlotFcns
            % function stop = outfun_d_opt(x,optimValues,state)
            %
            % switch state
            %     case 'iter'
            %         % Make updates to plot or guis as needed
            %         this_axis = gca;
            %         sensor_locations_k = vector_to_sensor_locations(x);
            %         img_k = assign_sensor_locations(imgi,sensor_locations_k);
            %         plot_sensors(img_k,false,'r','s',this_axis);
            %         axis([-1.1*R0 1.1*R0 -1.1*R0 1.1*R0 0 model_parameters.height])
            %         box on;grid on;
            %         view(3)
            %         drawnow
            %     case 'interrupt'
            %         % Probably no action here. Check conditions to see
            %         % whether optimization should quit.
            %     case 'init'
            %         hold on
            %         show_fem(imgi);
            %         % camlight;lighting gouraud
            %     case 'done'
            %         % Cleanup of plots, guis, or final plot
            %         this_axis = gca;
            %         sensor_locations_k = vector_to_sensor_locations(x);
            %         img_k = assign_sensor_locations(imgi,sensor_locations_k);
            %         plot_sensors(img_k,false,'b','s',this_axis);
            %         axis([-1.1*R0 1.1*R0 -1.1*R0 1.1*R0 0 model_parameters.height])
            %         box on;grid on;
            %         view(3)
            %         drawnow
            %     otherwise
            % end
            %
            % stop = false; %continue
            % end
            %
            % function stop = outfun_a_opt(x,optimValues,state)
            % switch state
            %     case 'iter'
            %
            %         % Make updates to plot or guis as needed
            %         this_axis = gca;
            %         sensor_locations_k = vector_to_sensor_locations(x);
            %         img_k = assign_sensor_locations(imgi,sensor_locations_k);
            %         plot_sensors(img_k,false,'r','o',this_axis);
            %         axis([-1.1*R0 1.1*R0 -1.1*R0 1.1*R0 0 model_parameters.height])
            %         box on;grid on;
            %         view(3)
            %         drawnow
            %     case 'interrupt'
            %         % Probably no action here. Check conditions to see
            %         % whether optimization should quit.
            %     case 'init'
            %         hold on
            %         show_fem(imgi);
            %         % camlight; lighting gouraud
            %
            %     case 'done'
            %         % Cleanup of plots, guis, or final plot
            %         this_axis = gca;
            %         sensor_locations_k = vector_to_sensor_locations(x);
            %         img_k = assign_sensor_locations(imgi,sensor_locations_k);
            %         plot_sensors(img_k,false,'b','x',this_axis);
            %         axis([-1.1*R0 1.1*R0 -1.1*R0 1.1*R0 0 model_parameters.height])
            %         box on;grid on;
            %         view(3)
            %         drawnow
            %     otherwise
            % end
            %
            % stop = false; %continue
            % end
            
            %% OutFcns
            function [xsol,fval,history] = runfmincon(funcgrad,q0,lb,ub,options)
                % Set up shared variables with outfun
                history.x = [];
                history.fval = [];
                
                options.OutputFcn = @outfun;

                function stop = outfun(x,optimValues,state)
                    stop = false;

                    switch state
                        case 'init'
                            history.fval = [history.fval, optimValues.fval];
                            history.x = [history.x, x];
                        case 'iter'
                            history.fval = [history.fval, optimValues.fval];
                            history.x = [history.x, x];
                        case 'done'
                        otherwise
                    end
                end

                [xsol,fval] = fmincon(funcgrad,q0,[],[],[],[],lb,ub,[],options);

            end
            %% Define function+gradient function

            function [func,grad] = funcwithgrad(q,f_impl,g_impl)
                % Calculate objective f
                func = f_impl(q);

                if nargout > 1 % gradient required
                    grad =  g_impl(q);
                end
            end
            

            %% Optimize 1-axis A-optimality with fmincon in z-coordinate of sensors
            fprintf('(Shell) 1-axis A-optimality OED - fmincon\n')
            
            options = optimoptions('fmincon',...
                'OptimalityTolerance',1e-9,'Display','iter','MaxIterations',max_iteratons, ...
                'Algorithm','interior-point','HessianApproximation',hessian_approximation,...
                'SpecifyObjectiveGradient',true,'UseParallel',use_parallel);

            fun_a_opt_mdeit_1 = @(q) funcwithgrad(q,f_a_opt_mdeit_dim_con,g_a_opt_mdeit_dim_con);

            lb = [ r0*ones(1,n_sensors),   [0,pi],   [0 0]];
            ub = [ r0*ones(1,n_sensors),   [0,pi],   model_parameters.height*ones(1,n_sensors) ];

            figure;        

            all_x_a_opt_con_shell = zeros(n_trial,length(x0));
            img_a_opt_con_shell = cell(1,n_trial);
            history = cell(1,n_trial);
            for n = 1:n_trial

                % Plot the initial sensor locations
                cla;
                show_fem(imgi);
                img_temp = assign_sensor_locations(imgi,all_sensor_locations_0(:,:,n));
                plot_sensors(img_temp,false,'b','.');
                hold off
                box on;grid on;
                axis([-1.1*rmax 1.1*rmax -1.1*rmax 1.1*rmax 0 model_parameters.height])
                drawnow

                % Initial conditions
                q0_con = sensor_locations_to_vector_con(all_sensor_locations_0(:,:,n));

                [x_a_opt_con_shell,fval,history{n}] = runfmincon(fun_a_opt_mdeit_1,q0_con,lb,ub,options);

                % x_a_opt_con_shell = fmincon(fun_a_opt_mdeit_1,q0_con,[],[],[],[],lb,ub,[],options);
                
                all_x_a_opt_con_shell(n,:) = x_a_opt_con_shell;
                img_a_opt_con_shell{n} = assign_sensor_locations(imgi,vector_to_sensor_locations_con(x_a_opt_con_shell));
            end

            %% Plot cost function over grid of degrees of freedom and the iterates
            z1 = linspace(0,model_parameters.height,10);
            z2 = z1;
            [Z1,Z2] = meshgrid(z1,z2);

            cost = zeros(size(Z1));

            q0 = [r0 r0 0 pi 0 0];
            for i = 1:size(Z1,1)
                fprintf('Working on line %i of %i\n',i,size(Z1,1));
                for j = 1:size(Z1,2)
                    q = q0;
                    q(5) = Z1(i,j);
                    q(6) = Z2(i,j);
                    cost(i,j) = f_a_opt_mdeit_dim_con(q);
                end
            end

            surf(Z1,Z2,cost);
            xlabel('$z_1$','interpreter','latex')
            ylabel('$z_1$','interpreter','latex')
            grid on;grid minor; box on;
            
            hold on
            colors = ['r','b'];
            for n = 1:n_trial
                n_iterates = numel(history{n}.fval);

                z1 = history{n}.x(5,1:n_iterates);
                z2 = history{n}.x(6,1:n_iterates);
                fz = history{n}.fval(1:n_iterates);

                text(z1,z2,fz,num2str((1:n_iterates)'),'Color',colors)
                plot3(z1,z2,fz,'.','MarkerSize',20,'Color',colors(n))
                drawnow
            end
            hold off

            %% Plot the optimal sensor locations
            theta = linspace(0,2*pi,100);
            [x, y, z] = cylinder(R0, 50);
            z = z * 3;


            figure
            hold on
            show_fem(imgi);
            
            for n = 1:n_trial
                img_init = assign_sensors_locations(imgi,all_sensor_locations_0(:,:,n));
                plot_sensors(img_init,false,'b','.');
            end

            for n = 1:n_trial
                img_temp = img_a_opt_con_shell{n};
                plot_sensors(img_temp,false,'r','s');
            end
            
            hold on
            hSurf = surf(h,x, y, z);
            set(hSurf, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor', 'red')
            
            plot3(h,R0*cos(theta),R0*sin(theta),model_parameters.height/2*ones(size(theta)),'LineStyle','--','Color','b')
            axis([-1.1*rmax 1.1*rmax -1.1*rmax 1.1*rmax 0 model_parameters.height])
            hold off
            box on;grid on;
            camlight; lighting gouraud
            drawnow
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

            model_parameters.maxsz = max(model_parameters.height,model_parameters.radius)/8;
            model_parameters.numOfElectrodesPerRing = 8;
            model_parameters.numOfRings = 4;
            model_parameters.numOfSensors = 4;
            model_parameters.sensorRadius = model_parameters.radius*1.5;
            model_parameters.material.type = 'spherical';
            model_parameters.material.position(3) = 3/4*model_parameters.height;

            background_conductivity = 3.28e-1/sigma0;
            anomaly_conductivity = background_conductivity/10;

            %% Simulation parameters

            do_3_axis = true;
            dim = 3;            %if doing 1-axis, which dimmension?

            %Optimizer parameters
            max_iteratons = 10;
            hessian_approximation = 'bfgs';
            use_parallel = true;
            algorithm = 'quasi-newton';

            rmax = 2*model_parameters.radius;
            rmin = model_parameters.radius*1.1;

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
            imgi = add_material_properties(imgh, [background_conductivity,anomaly_conductivity]);

            % show_fem(imgi);
            % plot_sensors(imgi);

            n_sensors = numel(imgi.fwd_model.sensors);
            n_stim = numel(imgi.fwd_model.stimulation);
            n_elem = size(imgi.fwd_model.elems,1);
            n_nodes = size(fmdl.nodes,1);

            A = @(x) M(imgi,x);

            %% Assume white noise and white prior

            B = fwd_solve_mdeit(imgi);
            max_B = max([abs(B.Bx(:));abs(B.By(:));abs(B.Bz(:))]);

            noise_std_deviation = max_B/10;

            % Set the noise variance with respect to the data magnitude
            variance_noise = noise_std_deviation^2;

            % Check if we are in a prior dominated regime, or if J(\sigma)
            % actually contributes to the gradient of the sensor positions

            % Gamma_post = (1/variance_noise J'*J + 1/variance_prior I
            % )^-1. The trace of Gamma_post is given by:

            % \sum_{i=1}^n_elem alpha/(1+\lambda_i\alpha/beta), where
            % \lambda_i are the eigenvalues of J'J, alpha is the variance
            % of the prior and beta is the variance of the noise

            % If most \lambda_i\alpha/beta << 1, we are on the prior
            % dominated regime, and the trace will be n_elem*\alpha.

            % If we are in a data dominated regime, many
            % \lambda_i\alpha/beta >> 1. In that case, we can approximate
            % the trace as:

            % (n_elem - r)\alpha + \sum_{i=1}^r \beta/\lambda_i, where r is
            % the number of eigenvalues that satisfy \lambda_i\alpha/beta
            % >> 1.

            % Since the eigenvalues of J'J decay rapidly, its hard to force
            % the data-oriented regime

            if exist('do_3_axis','var')
                if ~do_3_axis
                    Jdim = calc_jacobian_1axis_direct_fully_vectorized_local(imgi,A,dim);

                    coeff = 50;
                    variance_prior = coeff*variance_noise/eigs(Jdim'*Jdim,1);

                    % Check how many eigenvectors of J'J are in the data-dominated
                    % regime, as opposed to the prior-dominated regime
                    d = sum(eigs(Jdim'*Jdim,n_elem).*variance_prior/variance_noise>1);

                    fprintf('# lambda_i*alpha:beta>1 = %i (%i)\n',d,n_elem);

                    % White noise with zero mean and \mu variance
                    Gamma_noise = variance_noise.*speye(n_stim*n_sensors,n_stim*n_sensors);
                    Gamma_noise_inv = inv(Gamma_noise);

                    Gamma_prior = variance_prior.*speye(n_elem,n_elem);
                    Gamma_prior_inv = inv(Gamma_prior);

                else
                    Jx = calc_jacobian_1axis_direct_fully_vectorized_local(imgi,A,1);
                    Jy = calc_jacobian_1axis_direct_fully_vectorized_local(imgi,A,2);
                    Jz = calc_jacobian_1axis_direct_fully_vectorized_local(imgi,A,3);

                    J_3axis = [Jx;Jy;Jz];

                    coeff = 50;
                    variance_prior_3_axis = coeff*variance_noise/eigs(J_3axis'*J_3axis,1);

                    d_3axis = sum(eigs(J_3axis'*J_3axis,n_elem).*variance_prior_3_axis/variance_noise>1);

                    fprintf('# lambda_i*alpha:beta>1 = %i (%i)\n',d_3axis,n_elem);

                    Gamma_prior = variance_prior_3_axis.*speye(n_elem,n_elem);
                    Gamma_prior_inv = inv(Gamma_prior);

                    Gamma_noise = variance_noise.*speye(3*n_stim*n_sensors,3*n_stim*n_sensors);
                    Gamma_noise_inv = inv(Gamma_noise);
                end
            else
                error('Here');
            end

            alpha = 0; %controls the force of the repulsion!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            %% Perform optimization of sensor position

            % Initialize
            R0 = model_parameters.radius*1.5;
            theta_new = 0:2*pi/n_sensors:2*pi*(n_sensors-1)/n_sensors;
            r_new = R0*ones(1,n_sensors);
            z_new = model_parameters.height/2*ones(1,n_sensors);

            sensor_locations_0 = [(r_new.*cos(theta_new))',(r_new.*sin(theta_new))',z_new'];

            img_init = assign_sensor_locations(imgi,sensor_locations_0);

            plot_sensors(img_init);

            R0 = max(sqrt(sensor_locations_0(:,1).^2+sensor_locations_0(:,2).^2));
            z0 = max(sensor_locations_0(:,3));

            % These functions map sensor locations to sensor theta,eta

            % Reparametrize z so we can enforce it to be in [0,3]. Let
            % z = 3*sigmoid(nu), then z\in[0,3] and \nu is
            % unconstrained

            % Reparametrize r so we can enforce it to be in [rmin,rmax]. Let
            % r = rmin + (rmax-rmin)*sigmoid(xi), then r\in[rmin,rmax] and \xi is
            % unconstrained

            function x = sensor_locations_to_vector(sensor_locations)

                n_sensors = size(sensor_locations,1);

                x = zeros(3*n_sensors,1);

                for m = 1:n_sensors

                    % \xi
                    rm =  sqrt(sensor_locations(m,1)^2+sensor_locations(m,2)^2);
                    x(m) = -log((rmax-rmin)/(rm-rmin) -1);

                    % \theta
                    x(n_sensors + m) = atan2(sensor_locations(m,2),sensor_locations(m,1)); %thetam

                    % \eta
                    zm =  sensor_locations(m,3);
                    x(2*n_sensors+m) = log(zm/(3-zm)); %etam
                end

            end

            function sensor_locations = vector_to_sensor_locations(x)

                assert(mod(numel(x),3)==0);
                n_sensors = numel(x)/3;

                sensor_locations = zeros(n_sensors,3);

                for m = 1:n_sensors

                    xim = x(m);
                    thetam = x(m+n_sensors);
                    etam = x(m+2*n_sensors);

                    rm = rmin + (rmax-rmin)/(1+exp(-xim));
                    zm = 3/(1+exp(-etam));

                    sensor_locations(m,:) = [rm*cos(thetam),rm*sin(thetam),zm];
                end
            end

            % Sanity check
            s = vector_to_sensor_locations(sensor_locations_to_vector(sensor_locations_0));
            assert(norm(s-sensor_locations_0)<1e-5,'Unexpected');

            x0 = sensor_locations_to_vector(sensor_locations_0);

            if exist('do_3_axis','var')
                if ~do_3_axis
                    % Scale repulsion term to be of the same order as initial cost
                    f_a_opt_0 = compute_cost_function_a_opt(imgi,sensor_locations_0,Gamma_prior_inv,Gamma_noise_inv,A,dim);
                    f_d_opt_0 = compute_cost_function_d_opt(imgi,sensor_locations_0,Gamma_prior_inv,Gamma_noise_inv,A,dim);

                    G_0 = repulsion_cost(sensor_locations_0);

                end
            end

            % Sanity check
            if exist('do_3_axis','var')
                if ~do_3_axis
                    out_f_a_opt = f_a_opt(imgi,x0,Gamma_prior_inv,Gamma_noise_inv,A,dim);
                    disp(out_f_a_opt);
                end
            end

            % Compute cost function ( 1-axis MDEIT)
            function out = f_a_opt(img,x,inv_Gamma_prior,inv_Gamma_noise,A,dim)

                sensor_locations = vector_to_sensor_locations(x);

                phi_oed = compute_cost_function_a_opt(img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A,dim);

                G = repulsion_cost(sensor_locations);

                out = phi_oed+alpha*abs(f_a_opt_0/G_0)*G;
            end

            function out = f_d_opt(img,x,inv_Gamma_prior,inv_Gamma_noise,A,dim)

                sensor_locations = vector_to_sensor_locations(x);

                phi_oed= compute_cost_function_d_opt(img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A,dim);

                G = repulsion_cost(sensor_locations);

                out = phi_oed+alpha*abs(f_d_opt_0/G_0)*G;
            end

            % Cost function gradients (1-axis MDEIT). Allow movement in xi
            % theta and eta.
            function dphi = g_a_opt_optimized(img,x,Gamma_prior_inv,Gamma_noise_inv,A,dim)

                sigmoid = @(y) 1./(1+exp(-y));
                n_sensors = numel(img.fwd_model.sensors);

                x = x(:);
                xim = x(1:n_sensors);
                thetam = x(n_sensors+1:2*n_sensors);
                etam = x(2*n_sensors+1:3*n_sensors);

                rm = rmin + (rmax-rmin).*sigmoid(xim);
                zm = 3.*sigmoid(etam);

                sensor_locations = vector_to_sensor_locations(x);

                dphidp = compute_cost_function_gradient_a_opt_optimized(...
                    img,sensor_locations,Gamma_prior_inv,Gamma_noise_inv,A,dim);

                dphidx = dphidp(1,:);
                dphidy = dphidp(2,:);
                dphidz = dphidp(3,:);

                dphidr = cos(thetam)'.*dphidx+sin(thetam)'.*dphidy;
                dphidtheta = -rm'.*sin(thetam)'.*dphidx + rm'.*cos(thetam)'.*dphidy;

                drdxi = (rmax-rmin).*sigmoid(xim).'.*(1-sigmoid(xim)).';
                dzdeta = 3.*sigmoid(etam).'.*(1-sigmoid(etam)).';

                dphidxi = dphidr.*drdxi;
                dphideta = dphidz.*dzdeta;

                % % Derivatives of repulsion cost w.r.t cartesian coordinates
                % [dGx,dGy,dGz] = repulsion_grad_cartesian(sensor_locations);
                %
                % % Derivatives of repulsion w.r.t theta and eta
                % dGdtheta = -rm*sin(thetam).'.*dGx(:).' +rm*cos(thetam).'.*dGy(:).';
                % dGdeta = 3*dGz(:).'.*sigmoid(etam).'.*(1-sigmoid(etam)).';
                %
                % % Correction factor
                % dGdtheta = alpha*abs(f_a_opt_0/G_0)*dGdtheta;
                % dGdeta = alpha*abs(f_a_opt_0/G_0)*dGdeta;

                dphi = [dphidxi,dphidtheta,dphideta];
            end

            % For trying with fmincon
            %( 1-axis MDEIT)
            function out = f_a_opt_constrained(img,x,inv_Gamma_prior,inv_Gamma_noise,A,dim)

                assert(mod(numel(x),3)==0);
                n_sensors = numel(x)/3;

                sensor_locations = zeros(n_sensors,3);

                for m = 1:n_sensors
                    rm = x(m);
                    thetam = x(m+n_sensors);
                    zm = x(m+2*n_sensors);

                    sensor_locations(m,:) = [rm*cos(thetam),rm*sin(thetam),zm];
                end

                phi_oed = compute_cost_function_a_opt(img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A,dim);

                out = phi_oed;
            end
            function dphi = g_a_opt_optimized_constrained(img,x,Gamma_prior_inv,Gamma_noise_inv,A,dim)

                n_sensors = numel(img.fwd_model.sensors);

                x = x(:);
                rm = x(1:n_sensors);
                thetam = x(n_sensors+1:2*n_sensors);
                zm = x(2*n_sensors+1:3*n_sensors);

                sensor_locations = vector_to_sensor_locations(x);

                dphidp = compute_cost_function_gradient_a_opt_optimized(...
                    img,sensor_locations,Gamma_prior_inv,Gamma_noise_inv,A,dim);

                dphidx = dphidp(1,:);
                dphidy = dphidp(2,:);
                dphidz = dphidp(3,:);

                dphidr = cos(thetam)'.*dphidx+sin(thetam)'.*dphidy;
                dphidtheta = -rm'.*sin(thetam)'.*dphidx + rm'.*cos(thetam)'.*dphidy;

                dphi = [dphidr,dphidtheta,dphidz];
            end


            % ( 3-axis MDEIT)
            function out = f_a_opt_constrained_3_axis(img,x,inv_Gamma_prior,inv_Gamma_noise,A)

                assert(mod(numel(x),3)==0);
                n_sensors = numel(x)/3;

                sensor_locations = zeros(n_sensors,3);

                for m = 1:n_sensors
                    rm = x(m);
                    thetam = x(m+n_sensors);
                    zm = x(m+2*n_sensors);

                    sensor_locations(m,:) = [rm*cos(thetam),rm*sin(thetam),zm];
                end

                phi_oed = compute_cost_function_a_opt_3_axis(img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A);

                out = phi_oed;
            end
            function dphi = g_a_opt_optimized_constrained_3_axis(img,x,Gamma_prior_inv,Gamma_noise_inv,A)

                n_sensors = numel(img.fwd_model.sensors);

                x = x(:);
                rm = x(1:n_sensors);
                thetam = x(n_sensors+1:2*n_sensors);
                zm = x(2*n_sensors+1:3*n_sensors);

                sensor_locations = vector_to_sensor_locations(x);

                dphidp = compute_cost_function_gradient_a_opt_optimized_3_axis(...
                    img,sensor_locations,Gamma_prior_inv,Gamma_noise_inv,A);

                dphidx = dphidp(1,:);
                dphidy = dphidp(2,:);
                dphidz = dphidp(3,:);

                dphidr = cos(thetam)'.*dphidx+sin(thetam)'.*dphidy;
                dphidtheta = -rm'.*sin(thetam)'.*dphidx + rm'.*cos(thetam)'.*dphidy;

                dphi = [dphidr,dphidtheta,dphidz];
            end


            function dphi = g_d_opt_optimized(img,x,Gamma_prior_inv,Gamma_noise_inv,A,dim)

                sigmoid = @(y) 1./(1+exp(-y));
                n_sensors = numel(img.fwd_model.sensors);

                x = x(:);
                xim = x(1:n_sensors);
                thetam = x(n_sensors+1:2*n_sensors);
                etam = x(2*n_sensors+1:3*n_sensors);

                rm = rmin + (rmax-rmin).*sigmoid(xim);
                zm = 3.*sigmoid(etam);

                sensor_locations = vector_to_sensor_locations(x);

                % Derivatives of OED cost function w.r.t. cartesian
                % coordinates
                dphidp = compute_cost_function_gradient_d_opt_optimized(...
                    img,sensor_locations,Gamma_prior_inv,Gamma_noise_inv,A,dim);

                dphidx = dphidp(1,:);
                dphidy = dphidp(2,:);
                dphidz = dphidp(3,:);

                % Derivatives of OED cost function w.r.t. theta and eta
                dphidr = cos(thetam)'.*dphidx+sin(thetam)'.*dphidy;
                drdxi = (rmax-rmin).*sigmoid(xim).'.*(1-sigmoid(xim)).';
                dzdeta = 3.*sigmoid(etam).'.*(1-sigmoid(etam)).';

                dphidxi = dphidr.*drdxi;
                dphidtheta = -rm'.*sin(thetam)'.*dphidx + rm'.*cos(thetam)'.*dphidy;
                dphideta = dphidz.*dzdeta;

                % % Derivatives of repulsion cost w.r.t cartesian coordinates
                % [dGx,dGy,dGz] = repulsion_grad_cartesian(sensor_locations);
                %
                % % Derivatives of repulsion w.r.t theta and eta
                % dGdtheta = -rm*sin(thetam).'.*dGx(:).' +rm*cos(thetam).'.*dGy(:).';
                % dGdeta = 3*dGz(:).'.*sigmoid(etam).'.*(1-sigmoid(etam)).';
                %
                % % Correction factor
                % dGdtheta = alpha*abs(f_d_opt_0/G_0)*dGdtheta;
                % dGdeta = alpha*abs(f_d_opt_0/G_0)*dGdeta;

                dphi = [dphidxi,dphidtheta,dphideta];
            end

            % Compute cost function ( 3-axis MDEIT)
            function out = f_a_opt_3_axis(img,x,inv_Gamma_prior,inv_Gamma_noise,A)

                sensor_locations = vector_to_sensor_locations(x);

                phi_oed = compute_cost_function_a_opt_3_axis(img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A);

                out = phi_oed;
            end

            function out = f_d_opt_3_axis(img,x,inv_Gamma_prior,inv_Gamma_noise,A)

                sensor_locations = vector_to_sensor_locations(x);

                phi_oed= compute_cost_function_d_opt_3_axis(img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A);

                out = phi_oed;
            end

            % Cost function gradients (3-axis MDEIT). Allow movement in
            % theta and eta.
            function dphi = g_a_opt_optimized_3_axis(img,x,Gamma_prior_inv,Gamma_noise_inv,A)
                sigmoid = @(y) 1./(1+exp(-y));
                n_sensors = numel(img.fwd_model.sensors);

                x = x(:);
                rm = R0;
                thetam = x(1:n_sensors);
                etam = x(n_sensors+1:2*n_sensors);

                sensor_locations = vector_to_sensor_locations(x);

                dphidp = compute_cost_function_gradient_a_opt_optimized_3_axis(...
                    img,sensor_locations,Gamma_prior_inv,Gamma_noise_inv,A);

                dphidx = dphidp(1,:);
                dphidy = dphidp(2,:);
                dphidz = dphidp(3,:);

                dphidtheta = -rm'.*sin(thetam)'.*dphidx + rm'.*cos(thetam)'.*dphidy;
                dphideta = 3*dphidz.*sigmoid(etam).'.*(1-sigmoid(etam)).';

                dphi = [dphidtheta,dphideta];
            end

            function dphi = g_d_opt_optimized_3_axis(img,x,Gamma_prior_inv,Gamma_noise_inv,A)
                sigmoid = @(y) 1./(1+exp(-y));
                n_sensors = numel(img.fwd_model.sensors);

                x = x(:);
                rm = R0;
                thetam = x(1:n_sensors);
                etam = x(n_sensors+1:2*n_sensors);
                sensor_locations = vector_to_sensor_locations(x);

                dphidp = compute_cost_function_gradient_d_opt_optimized_3_axis(...
                    img,sensor_locations,Gamma_prior_inv,Gamma_noise_inv,A);

                dphidx = dphidp(1,:);
                dphidy = dphidp(2,:);
                dphidz = dphidp(3,:);

                dphidtheta = -rm'.*sin(thetam)'.*dphidx + rm'.*cos(thetam)'.*dphidy;
                dphideta = 3*dphidz.*sigmoid(etam).'.*(1-sigmoid(etam)).';

                dphi = [dphidtheta,dphideta];
            end

            %% PlotFcns

            function stop = outfun_d_opt(x,optimValues,state)
                switch state
                    case 'iter'
                        % Make updates to plot or guis as needed
                        this_axis = gca;
                        sensor_locations_k = vector_to_sensor_locations(x);
                        img_k = assign_sensor_locations(imgi,sensor_locations_k);
                        plot_sensors(img_k,false,'r','s',this_axis);
                        axis([-1.1*R0 1.1*R0 -1.1*R0 1.1*R0 0 model_parameters.height])
                        box on;grid on;
                        view(3)
                        drawnow
                    case 'interrupt'
                        % Probably no action here. Check conditions to see
                        % whether optimization should quit.
                    case 'init'
                        hold on
                        show_fem(imgi);
                        % camlight;lighting gouraud
                    case 'done'
                        % Cleanup of plots, guis, or final plot
                        this_axis = gca;
                        sensor_locations_k = vector_to_sensor_locations(x);
                        img_k = assign_sensor_locations(imgi,sensor_locations_k);
                        plot_sensors(img_k,false,'b','s',this_axis);
                        axis([-1.1*R0 1.1*R0 -1.1*R0 1.1*R0 0 model_parameters.height])
                        box on;grid on;
                        view(3)
                        drawnow
                    otherwise
                end

                stop = false; %continue
            end

            function stop = outfun_a_opt(x,optimValues,state)
                switch state
                    case 'iter'

                        % Make updates to plot or guis as needed
                        this_axis = gca;
                        sensor_locations_k = vector_to_sensor_locations(x);
                        img_k = assign_sensor_locations(imgi,sensor_locations_k);
                        plot_sensors(img_k,false,'r','o',this_axis);
                        axis([-1.1*R0 1.1*R0 -1.1*R0 1.1*R0 0 model_parameters.height])
                        box on;grid on;
                        view(3)
                        drawnow
                    case 'interrupt'
                        % Probably no action here. Check conditions to see
                        % whether optimization should quit.
                    case 'init'
                        hold on
                        show_fem(imgi);
                        % camlight; lighting gouraud

                    case 'done'
                        % Cleanup of plots, guis, or final plot
                        this_axis = gca;
                        sensor_locations_k = vector_to_sensor_locations(x);
                        img_k = assign_sensor_locations(imgi,sensor_locations_k);
                        plot_sensors(img_k,false,'b','x',this_axis);
                        axis([-1.1*R0 1.1*R0 -1.1*R0 1.1*R0 0 model_parameters.height])
                        box on;grid on;
                        view(3)
                        drawnow
                    otherwise
                end

                stop = false; %continue
            end

            %% Use fmincon

            x0_con = zeros(3*n_sensors,1);

            for m = 1:n_sensors
                rm =  sqrt(sensor_locations_0(m,1)^2+sensor_locations_0(m,2)^2);
                thetam = atan2(sensor_locations_0(m,2),sensor_locations_0(m,1));
                zm = sensor_locations_0(m,3);

                x0_con(m) = rm;
                x0_con(m+n_sensors) = thetam;
                x0_con(m+2*n_sensors) = zm;
            end

            options = optimoptions('fmincon',...
                'OptimalityTolerance',1e-9,'Display','iter','MaxIterations',max_iteratons, ...
                'Algorithm','interior-point','HessianApproximation',hessian_approximation,...
                'SpecifyObjectiveGradient',true,'UseParallel',use_parallel,...
                PlotFcn=@outfun_a_opt);

            if exist('do_3_axis','var')
                if do_3_axis
                    fprintf('Doing 3-axis A-optimality OED - fmincon\n')
                    fun_a_opt_constrained = @funcwithgrad_a_opt_constrained_3_axis;
                else
                    fprintf('Doing 1-axis A-optimality OED - fmincon\n')
                    fun_a_opt_constrained = @funcwithgrad_a_opt_constrained;
                end
            else
                error('Here');
            end

            function [f,g] = funcwithgrad_a_opt_constrained(x)
                % Calculate objective f
                f = f_a_opt_constrained(imgi,x,Gamma_prior_inv,Gamma_noise_inv,A,dim);

                if nargout > 1 % gradient required
                    g =  g_a_opt_optimized_constrained(imgi,x,Gamma_prior_inv,Gamma_noise_inv,A,dim);
                end
            end

            function [f,g] = funcwithgrad_a_opt_constrained_3_axis(x)
                % Calculate objective f
                f = f_a_opt_constrained_3_axis(imgi,x,Gamma_prior_inv,Gamma_noise_inv,A);

                if nargout > 1 % gradient required
                    g =  g_a_opt_optimized_constrained_3_axis(imgi,x,Gamma_prior_inv,Gamma_noise_inv,A);
                end
            end

            lb = [ rmin*ones(1,n_sensors),   -inf(1,n_sensors),   0*ones(1,n_sensors)];
            ub = [ rmax*ones(1,n_sensors),    inf(1,n_sensors),   model_parameters.height*ones(1,n_sensors) ];

            x_a_opt_con = fmincon(fun_a_opt_constrained,x0_con,[],[],[],[],lb,ub,[],options);

            sensor_locations_a_opt_con = zeros(n_sensors,3);

            for m = 1:n_sensors
                rm = x_a_opt_con(m);
                thetam = x_a_opt_con(m+n_sensors);
                zm = x_a_opt_con(m+2*n_sensors);

                sensor_locations_a_opt_con(m,:) = [rm*cos(thetam),rm*sin(thetam),zm];
            end

            img_a_opt_con = assign_sensor_locations(imgi,sensor_locations_a_opt_con);

            %% Optimize A-optimality with trust region method native to MATLAB
            options = optimoptions('fminunc',...
                'OptimalityTolerance',1e-9,'Display','iter','MaxIterations',max_iteratons, ...
                'Algorithm', algorithm,'HessianApproximation',hessian_approximation,...
                'SpecifyObjectiveGradient',true,'UseParallel',use_parallel,...
                PlotFcn=@outfun_a_opt);

            if exist('do_3_axis','var')
                if do_3_axis
                    fprintf('Doing 3-axis A-optimality OED\n')
                    fun_a_opt = @funcwithgrad_a_opt_3_axis;
                else
                    fprintf('Doing 1-axis A-optimality OED\n')
                    fun_a_opt = @funcwithgrad_a_opt;
                end
            else
                error('Here');
            end

            function [f,g] = funcwithgrad_a_opt(x)
                % Calculate objective f
                f = f_a_opt(imgi,x,Gamma_prior_inv,Gamma_noise_inv,A,dim);

                if nargout > 1 % gradient required
                    g =  g_a_opt_optimized(imgi,x,Gamma_prior_inv,Gamma_noise_inv,A,dim);
                end
            end

            function [f,g] = funcwithgrad_a_opt_3_axis(x)
                % Calculate objective f
                f = f_a_opt_3_axis(imgi,x,Gamma_prior_inv,Gamma_noise_inv,A);

                if nargout > 1 % gradient required
                    g =  g_a_opt_optimized_3_axis(imgi,x,Gamma_prior_inv,Gamma_noise_inv,A);
                end
            end


            x_a_opt = fminunc(fun_a_opt,x0,options);
            sensor_locations_a_opt = vector_to_sensor_locations(x_a_opt);
            img_a_opt = assign_sensor_locations(imgi,sensor_locations_a_opt);

            % Disp proportion of costs
            if exist('do_3_axis','var')
                if ~do_3_axis
                    phi_oed_a_0 = f_a_opt(imgi, x0 ,Gamma_prior_inv,Gamma_noise_inv,A,dim);
                    phi_oed_a_opt= f_a_opt(imgi, x_a_opt ,Gamma_prior_inv,Gamma_noise_inv,A,dim);

                    G_a_opt = repulsion_cost(sensor_locations_a_opt);
                    total_cost_a_opt = phi_oed_a_opt + alpha*abs(f_a_opt_0/G_0)*G_a_opt;

                    fprintf('Total cost: %2.2g\n',total_cost_a_opt);
                    fprintf('OED cost: %2.2g -> %2.2g (%1.2f %%)\n',...se
                        phi_oed_a_0,phi_oed_a_opt,...
                        sign(phi_oed_a_opt-phi_oed_a_0)*abs(phi_oed_a_opt-phi_oed_a_0)/abs(phi_oed_a_0)*100);
                    fprintf('Repulsion cost: %2.2g -> %2.2g (%1.2f %%)\n',...
                        G_0,G_a_opt,...
                        sign(G_a_opt-G_0)*abs(G_a_opt-G_0)/abs(G_0)*100);
                else
                    phi_oed_a_0 = f_a_opt_3_axis(imgi, x0 ,Gamma_prior_inv,Gamma_noise_inv,A);
                    phi_oed_a_opt= f_a_opt_3_axis(imgi, x_a_opt ,Gamma_prior_inv,Gamma_noise_inv,A);

                    total_cost_a_opt = phi_oed_a_opt;

                    fprintf('Total cost: %2.2g\n',total_cost_a_opt);
                    fprintf('OED cost: %2.2g -> %2.2g (%1.2f %%)\n',...se
                        phi_oed_a_0,phi_oed_a_opt,...
                        sign(phi_oed_a_opt-phi_oed_a_0)*abs(phi_oed_a_opt-phi_oed_a_0)/abs(phi_oed_a_0)*100);
                end
            end

            %% Optimize D-optimality with trust region method native to MATLAB

            options = optimoptions('fminunc',...
                'OptimalityTolerance',1e-3,'Display','iter','MaxIterations',max_iteratons,...
                'Algorithm', algorithm,'HessianApproximation',hessian_approximation,...
                'SpecifyObjectiveGradient',true,'UseParallel',use_parallel,...
                PlotFcn=@outfun_d_opt);
            % options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton','SpecifyObjectiveGmax_iteratonsradient',true);

            if exist('do_3_axis','var')
                if do_3_axis
                    fprintf('Doing 3-axis D-optimality OED\n')
                    fun_d_opt = @funcwithgrad_d_opt_3_axis;
                else
                    fprintf('Doing 1-axis D-optimality OED\n')
                    fun_d_opt = @funcwithgrad_d_opt;
                end
            else
                error('Here');
            end

            function [f,g] = funcwithgrad_d_opt(x)
                % Calculate objective f
                f = f_d_opt(imgi,x,Gamma_prior_inv,Gamma_noise_inv,A,dim);

                if nargout > 1 % gradient required
                    g =  g_d_opt_optimized(imgi,x,Gamma_prior_inv,Gamma_noise_inv,A,dim);
                end
            end

            function [f,g] = funcwithgrad_d_opt_3_axis(x)
                % Calculate objective f
                f = f_d_opt_3_axis(imgi,x,Gamma_prior_inv,Gamma_noise_inv,A);

                if nargout > 1 % gradient required
                    g =  g_d_opt_optimized_3_axis(imgi,x,Gamma_prior_inv,Gamma_noise_inv,A);
                end
            end

            x_d_opt = fminunc(fun_d_opt,x0,options);
            sensor_locations_d_opt = vector_to_sensor_locations(x_d_opt);
            img_d_opt = assign_sensor_locations(imgi,sensor_locations_d_opt);

            % Disp proportion of costs
            if exist('do_3_axis','var')
                if ~do_3_axis
                    phi_oed_d_0 = f_d_opt(imgi,x0,Gamma_prior_inv,Gamma_noise_inv,A,dim);
                    phi_oed_d_opt= f_d_opt(imgi,x_d_opt,Gamma_prior_inv,Gamma_noise_inv,A,dim);

                    G_d_opt = repulsion_cost(sensor_locations_d_opt);
                    total_cost_d_opt = phi_oed_d_opt+alpha*abs(f_d_opt_0/G_0)*G_d_opt;

                    fprintf('Total cost: %2.2g\n',total_cost_d_opt);
                    fprintf('OED cost: %2.2g -> %2.2g (%1.2f %%)\n',...
                        phi_oed_d_0,phi_oed_d_opt,...
                        sign(phi_oed_d_opt-phi_oed_d_0)*abs(phi_oed_d_opt-phi_oed_d_0)/abs(phi_oed_d_0)*100);
                    fprintf('Repulsion cost: %2.2g -> %2.2g (%1.2f %%)\n',...
                        G_0,G_d_opt,...
                        sign(G_d_opt-G_0)*abs(G_d_opt-G_0)/abs(G_0)*100);
                else
                    phi_oed_d_0 = f_d_opt_3_axis(imgi,x0,Gamma_prior_inv,Gamma_noise_inv,A);
                    phi_oed_d_opt= f_d_opt_3_axis(imgi,x_d_opt,Gamma_prior_inv,Gamma_noise_inv,A);

                    total_cost_d_opt = phi_oed_d_opt;

                    fprintf('Total cost: %2.2g\n',total_cost_d_opt);
                    fprintf('OED cost: %2.2g -> %2.2g (%1.2f %%)\n',...
                        phi_oed_d_0,phi_oed_d_opt,...
                        sign(phi_oed_d_opt-phi_oed_d_0)*abs(phi_oed_d_opt-phi_oed_d_0)/abs(phi_oed_d_0)*100);
                end
            end


            %% Plot the optimal sensor locations
            theta = linspace(0,2*pi,100);
            [x, y, z] = cylinder(rmin, 50);
            [x2, y2, z2] = cylinder(rmax, 50);
            z = z * 3;
            z2 = z2*3;


            figure
            hold on
            show_fem(imgi);
            plot_sensors(img_init);
            plot_sensors(img_a_opt,false,'g','s');
            plot_sensors(img_d_opt,false,'r','o');
            h = plot_sensors(img_a_opt_con,false,'k','s');

            hold on
            hSurf = surf(h,x, y, z);
            set(hSurf, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor', 'red')
            hSurf2 = surf(h,x2, y2, z2);
            set(hSurf2, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor', 'blue')
            plot3(h,R0*cos(theta),R0*sin(theta),model_parameters.height/2*ones(size(theta)),'LineStyle','--','Color','b')
            axis([-1.1*rmax 1.1*rmax -1.1*rmax 1.1*rmax 0 model_parameters.height])
            hold off
            box on;grid on;
            camlight; lighting gouraud
            drawnow
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

function dR = dRmkj_xyz(fmdl,j)
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

% --- Compute dRdp for all sensors, and for all p ---
dRdx = zeros(numSensors,numElements);
dRdy = zeros(numSensors,numElements);
dRdz = zeros(numSensors,numElements);



for m = 1:numSensors
    rm = fmdl.sensors(m).position;  % 1 x 3

    % Evaluate quadrature points in all elements
    % xi: numQuadPts x 3 x numElements
    xi = reshape(v1', [3,1,numElements]) + ...
        J(:,1,:).*reshape(coord(:,1),[1,numQuadPts,1]) + ...
        J(:,2,:).*reshape(coord(:,2),[1,numQuadPts,1]) + ...
        J(:,3,:).*reshape(coord(:,3),[1,numQuadPts,1]); % 3 x numQuadPts x numElements

    xi = permute(xi,[2 1 3]); % numQuadPts x 3 x numElements

    % dm_vec = rm - xi, same size
    dm_vec = rm - xi; % numQuadPts x 3 x numElements

    dm_norm2 = sum(dm_vec.^2,2);           % numQuadPts x 1 x numElements
    dm_j    = dm_vec(:,j,:);               % numQuadPts x 1 x numElements

    for p = 1:3
        % Compute fun values at all quadrature points
        % ((j==p)*norm(dm)^2 - 3*dm_j*dm_p) / norm(dm)^5

        dm_p    = dm_vec(:,p,:);

        funvals = ((j==p)*dm_norm2 - 3*dm_j.*dm_p) ./ (dm_norm2.^(5/2)); % numQuadPts x 1 x numElements

        % Integrate over quadrature points using weights
        switch p
            case 1
                dRdx(m,:) = squeeze(sum(funvals .* reshape(weights, [numQuadPts,1,1]),1))' .* (detJ/6);
            case 2
                dRdy(m,:) = squeeze(sum(funvals .* reshape(weights, [numQuadPts,1,1]),1))' .* (detJ/6);
            case 3
                dRdz(m,:) = squeeze(sum(funvals .* reshape(weights, [numQuadPts,1,1]),1))' .* (detJ/6);
        end
    end
end

dR = {dRdx,dRdy,dRdz};

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


function [dlambda,dR1,dR2] = compute_dlambda_xyz(img,dim)

num_sensors = numel(img.fwd_model.sensors);
n_nodes = size(img.fwd_model.nodes,1);

dlambdaX = zeros(n_nodes,num_sensors);
dlambdaY = zeros(n_nodes,num_sensors);
dlambdaZ = zeros(n_nodes,num_sensors);

mu_factor = img.fwd_model.mu0/(4*pi);

G = img.fwd_model.G;

switch dim
    case 1
        dR1 = dRmkj_xyz(img.fwd_model, 3);
        dR2 = dRmkj_xyz(img.fwd_model, 2);

        G1 = G.Gy;
        G2 = G.Gz;
    case 2
        dR1  = dRmkj_xyz(img.fwd_model, 1);
        dR2 = dRmkj_xyz(img.fwd_model, 3);

        G1 = G.Gz;
        G2 = G.Gx;
    case 3
        dR1 = dRmkj_xyz(img.fwd_model, 2);
        dR2 = dRmkj_xyz(img.fwd_model, 1);

        G1 = G.Gx;
        G2 = G.Gy;
    otherwise
        error('Invalid dimension');
end

sigma = img.elem_data;
A_matrix = M(img,sigma);

% Jacobi preconditioner - matrix free
d = sqrt(diag(A_matrix));        % vector of diagonal entries

Mfun = @(x) x ./ d;              % left preconditioner  M^{-1} x
Nfun = @(x) x ./ d;              % right preconditioner

% Multiply elementwise by sigma (diagonal) without creating sparse diagonal matrix
% rhs = mu_factor*( dR1 * Sigma * G1 -  dR2* Sigma * G2);
rhs1 = mu_factor*( (dR1{1} .* sigma.')*G1 - (dR2{1} .* sigma.')*G2 );
rhs2 = mu_factor*( (dR1{2} .* sigma.')*G1 - (dR2{2} .* sigma.')*G2 );
rhs3 = mu_factor*( (dR1{3} .* sigma.')*G1 - (dR2{3} .* sigma.')*G2 );

parfor m = 1:num_sensors

    [dlambdaX(:,m),~,~] = pcg(A_matrix,rhs1(m,:)',1e-9,numel(sigma),Mfun,Nfun);
    [dlambdaY(:,m),~,~] = pcg(A_matrix,rhs2(m,:)',1e-9,numel(sigma),Mfun,Nfun);
    [dlambdaZ(:,m),~,~] = pcg(A_matrix,rhs3(m,:)',1e-9,numel(sigma),Mfun,Nfun);
end

dlambda = {dlambdaX,dlambdaY,dlambdaZ};

end


function [dlambda,dR1t,dR2t] = compute_dlambdaxyz_xyz(img)

num_sensors = numel(img.fwd_model.sensors);
n_nodes = size(img.fwd_model.nodes,1);

dlambda = cell(3,3);

dlambdaX = zeros(n_nodes,num_sensors);
dlambdaY = zeros(n_nodes,num_sensors);
dlambdaZ = zeros(n_nodes,num_sensors);

dR1t = cell(3,3);
dR2t = cell(3,3);


mu_factor = img.fwd_model.mu0/(4*pi);

G = img.fwd_model.G;

sigma = img.elem_data;
A_matrix = M(img,sigma);

% Jacobi preconditioner - matrix free
d = sqrt(diag(A_matrix));        % vector of diagonal entries

Mfun = @(x) x ./ d;              % left preconditioner  M^{-1} x
Nfun = @(x) x ./ d;              % right preconditioner

for dim = 1:3

    switch dim
        case 1
            dR1 = dRmkj_xyz(img.fwd_model, 3);
            dR2 = dRmkj_xyz(img.fwd_model, 2);

            G1 = G.Gy;
            G2 = G.Gz;
        case 2
            dR1  = dRmkj_xyz(img.fwd_model, 1);
            dR2 = dRmkj_xyz(img.fwd_model, 3);

            G1 = G.Gz;
            G2 = G.Gx;
        case 3
            dR1 = dRmkj_xyz(img.fwd_model, 2);
            dR2 = dRmkj_xyz(img.fwd_model, 1);

            G1 = G.Gx;
            G2 = G.Gy;
    end

    % Multiply elementwise by sigma (diagonal) without creating sparse diagonal matrix
    % rhs = mu_factor*( dR1 * Sigma * G1 -  dR2* Sigma * G2);

    rhs1 = mu_factor*( (dR1{1} .* sigma.')*G1 - (dR2{1} .* sigma.')*G2 ); %for p=1
    rhs2 = mu_factor*( (dR1{2} .* sigma.')*G1 - (dR2{2} .* sigma.')*G2 ); %for p=2
    rhs3 = mu_factor*( (dR1{3} .* sigma.')*G1 - (dR2{3} .* sigma.')*G2 ); %for p=3

    parfor m = 1:num_sensors

        [dlambdaX(:,m),~,~] = pcg(A_matrix,rhs1(m,:)',1e-9,numel(sigma),Mfun,Nfun);
        [dlambdaY(:,m),~,~] = pcg(A_matrix,rhs2(m,:)',1e-9,numel(sigma),Mfun,Nfun);
        [dlambdaZ(:,m),~,~] = pcg(A_matrix,rhs3(m,:)',1e-9,numel(sigma),Mfun,Nfun);
    end

    dlambda{dim,1} = dlambdaX;
    dlambda{dim,2} = dlambdaY;
    dlambda{dim,3} = dlambdaZ;

    dR1t{dim,1} = dR1{1};
    dR1t{dim,2} = dR1{2};
    dR1t{dim,3} = dR1{3};

    dR2t{dim,1} = dR2{1};
    dR2t{dim,2} = dR2{2};
    dR2t{dim,3} = dR2{3};
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

function [dJ,dJ1,dJ2] = compute_dJ_xyz(img,dim)
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
[dlambda,dR1dp,dR2dp] = compute_dlambda_xyz(img,dim);

% Compute derivatives of R w.r.t sensor positions
switch dim
    case 1
        G1U = GyU;
        G2U = GzU;
    case 2
        G1U = GzU;
        G2U = GxU;
    case 3
        G1U = GxU;
        G2U = GyU;
    otherwise
        error('Dimension not supported.');
end

% Preallocate
dJ = cell(1,3);
dJ1 = zeros(num_sensors,num_stim,n_elem);
dJ2 = zeros(num_sensors,num_stim,n_elem);

%% Loop only over sensors
for p = 1:3
    for m = 1:num_sensors
        %% --- dJ1: contribution from dlambda ---
        % dlambda * G matrices, size: n_elem x 1
        dlGx = G.Gx*dlambda{p}(:,m);
        dlGy = G.Gy*dlambda{p}(:,m);
        dlGz = G.Gz*dlambda{p}(:,m);

        % Elementwise multiplication and sum over components
        % tmp: n_elem x num_stim
        tmp = dlGx .* GxU + dlGy .* GyU + dlGz .* GzU;

        % Multiply by element volumes, permute to [num_stim x n_elem]
        dJ1(m,:,:) = tmp.' .* elemV(:).';

        %% --- dJ2: contribution from dR/dp ---
        % dRydp, dRzdp: 1 x n_elem
        % Multiply by G*u per stimulation, permute to [num_stim x n_elem]
        dJ2(m,:,:) = -mu_factor * ( ...
            dR2dp{p}(m,:) .* G2U.' - dR1dp{p}(m,:) .* G1U.' ...
            );
    end

    %% Total derivative
    dJ{p} = dJ1 - dJ2;
end

end

function dJ = compute_dJxyz_xyz(img)
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
[dlambda,dR1dp,dR2dp] = compute_dlambdaxyz_xyz(img);

% Compute derivatives of R w.r.t sensor positions
G1U = cell(1,3);
G2U = cell(1,3);

for dim = 1:3
    switch dim
        case 1
            G1U{dim} = GyU;
            G2U{dim} = GzU;
        case 2
            G1U{dim} = GzU;
            G2U{dim} = GxU;
        case 3
            G1U{dim} = GxU;
            G2U{dim} = GyU;
        otherwise
            error('Dimension not supported.');
    end
end

% Preallocate
dJ = cell(3,3);

dJ1 = zeros(num_sensors,num_stim,n_elem);
dJ2 = zeros(num_sensors,num_stim,n_elem);


%% Loop only over sensors
for dim = 1:3
    for p = 1:3
        for m = 1:num_sensors
            %% --- dJ1: contribution from dlambda ---
            % dlambda * G matrices, size: n_elem x 1
            dlGx = G.Gx*dlambda{p}(:,m);
            dlGy = G.Gy*dlambda{p}(:,m);
            dlGz = G.Gz*dlambda{p}(:,m);

            % Elementwise multiplication and sum over components
            % tmp: n_elem x num_stim
            tmp = dlGx .* GxU + dlGy .* GyU + dlGz .* GzU;

            % Multiply by element volumes, permute to [num_stim x n_elem]
            dJ1(m,:,:) = tmp.' .* elemV(:).';

            %% --- dJ2: contribution from dR/dp ---
            switch dim
                case 1
                    dJ2(m,:,:) = -mu_factor * ( ...
                        dR2dp{dim,p}(m,:) .* G2U{dim}.' - dR1dp{dim,p}(m,:) .* G1U{dim}.' ...
                        );
                case 2
                    dJ2(m,:,:) = -mu_factor * ( ...
                        dR2dp{dim,p}(m,:) .* G2U{dim}.' - dR1dp{dim,p}(m,:) .* G1U{dim}.' ...
                        );
                case 3
                    dJ2(m,:,:) = -mu_factor * ( ...
                        dR2dp{dim,p}(m,:) .* G2U{dim}.' - dR1dp{dim,p}(m,:) .* G1U{dim}.' ...
                        );
            end


        end

        %% Total derivative
        dJ{dim,p} = dJ1 - dJ2;
    end
end
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

% For 1 axis-MDEIT
function cost = compute_cost_function_d_opt(img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A,dim)

n_sensors = numel(img.fwd_model.sensors);
n_stim = numel(img.fwd_model.stimulation);
n_elem = size(img.fwd_model.elems,1);

% Assign sensor locations
img = assign_sensor_locations(img,sensor_locations);

% Compute the jacobian at current sensor locations
J = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,dim);

% Define the inverse posterior covariance matrix
H = J.'*inv_Gamma_noise*J+inv_Gamma_prior;

L = chol(H,'lower');
logdetH = 2*sum(log(diag(L)));

cost = -logdetH;
end

function cost = compute_cost_function_a_opt(img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A,dim)

n_sensors = numel(img.fwd_model.sensors);
n_stim = numel(img.fwd_model.stimulation);
n_elem = size(img.fwd_model.elems,1);

% Assign sensor locations
img = assign_sensor_locations(img,sensor_locations);

% Compute the jacobian at current sensor locations
J = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,dim);

% Define the inverse posterior covariance matrix
H = J.'*inv_Gamma_noise*J+inv_Gamma_prior;

% L = chol(H,'lower');
% cost = sum(1./diag(L).^2);

cost = trace(inv(H));
end


% For 3 axis-MDEIT
function cost = compute_cost_function_d_opt_3_axis(img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A)
% Assign sensor locations
img = assign_sensor_locations(img,sensor_locations);

% Compute the jacobian at current sensor locations
Jx = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,1);
Jy = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,2);
Jz = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,3);

J = [Jx;Jy;Jz];

% Define the inverse posterior covariance matrix
H = J.'*inv_Gamma_noise*J+inv_Gamma_prior;

L = chol(H,'lower');
logdetH = 2*sum(log(diag(L)));

cost = -logdetH;
end

function cost = compute_cost_function_a_opt_3_axis(img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A)

% Assign sensor locations
img = assign_sensor_locations(img,sensor_locations);

% Compute the jacobian at current sensor locations
Jx = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,1);
Jy = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,2);
Jz = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,3);

J = [Jx;Jy;Jz];

% Define the inverse posterior covariance matrix
H = J.'*inv_Gamma_noise*J+inv_Gamma_prior;

% L = chol(H,'lower');
% cost = sum(1./diag(L).^2);

cost = trace(inv(H));
end





function G = repulsion_cost(sensor_locations)

n_sensors = size(sensor_locations,1);

% --- pairwise differences ---
% Expand coordinates
X = sensor_locations(:,1);
Y = sensor_locations(:,2);
Z = sensor_locations(:,3);

DX = X - X.';    % [n_sensors x n_sensors]
DY = Y - Y.';
DZ = Z - Z.';

% Squared distances
D2 = DX.^2 + DY.^2 + DZ.^2+ 1e-12; %smooth correction

% --- Cost ---
invD2 = 1 ./ D2;

% remove diagonal
invD2(1:n_sensors+1:end) = 0;

% sum only i<j. We sum repeated pairs (i,j) and (j,i), so we multiply by
% 0.5
G = 0.5 * sum(invD2(:));

end

function [dGx,dGy,dGz] = repulsion_grad_cartesian(sensor_locations)

n_sensors = size(sensor_locations,1);

% --- pairwise differences ---
% Expand coordinates
X = sensor_locations(:,1);
Y = sensor_locations(:,2);
Z = sensor_locations(:,3);

DX = X - X.';    % [n_sensors x n_sensors]
DY = Y - Y.';
DZ = Z - Z.';

% Squared distances
D2 = DX.^2 + DY.^2 + DZ.^2+ 1e-12; %smooth correction

% --- Cost ---
invD2 = 1 ./ D2;

% remove diagonal
invD2(1:n_sensors+1:end) = 0;

% We need 1/d^4
invD4 = invD2.^2;

% Compute force matrix coefficients
C = -2 * invD4;

% zero diagonal
C(1:n_sensors+1:end) = 0;

% Each gradient component is:
% dG_i = sum_j C_ij * (p_i - p_j)

dGx = sum(C .* DX, 2);
dGy = sum(C .* DY, 2);
dGz = sum(C .* DZ, 2);

end



function dphidp = compute_cost_function_gradient_d_opt(...
    img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A,dim,p)

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
H = J.'*inv_Gamma_noise*J+inv_Gamma_prior;

% Intermediate variable W (different in a-opt to d-opt)
W =  inv_Gamma_noise*J*(inv(H)); %[n_stim*n_sensors,n_elem]

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

function dphidp = compute_cost_function_gradient_a_opt(...
    img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A,dim,p)

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
H = J.'*inv_Gamma_noise*J+inv_Gamma_prior;

% Intermediate variable W (different in a-opt to d-opt)
W =  inv_Gamma_noise*J*(inv(H)^2); %[n_stim*n_sensors,n_elem]

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


% For 1-axis MDEIT
function dphidp = compute_cost_function_gradient_d_opt_optimized(...
    img,sensor_locations,Gamma_prior_inv,Gamma_noise_inv,A,dim)

n_sensors = numel(img.fwd_model.sensors);
n_stim = numel(img.fwd_model.stimulation);
n_elem = size(img.fwd_model.elems,1);

img = assign_sensor_locations(img,sensor_locations);

J = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,dim);
dJds = compute_dJ_xyz(img,dim);

H = J.'*Gamma_noise_inv*J + Gamma_prior_inv;

L = chol(H,'lower');

% Compute H^{-2} through Cholesky factorization and linear system solves

% Solve L'Y = I -> LH^{-1} = Y -> L'W = H^{-1} -> L H^{-2} = W
Y = L'\eye(n_elem);
invH = L\Y;

W = Gamma_noise_inv*J*invH;

dphidp = zeros(3,n_sensors);

for p = 1:3
    for m = 1:n_sensors
        ids = m : n_sensors : n_sensors*n_stim;

        dJm = reshape(dJds{p}(m,:,:), [n_stim, n_elem]);
        Wm  = W(ids,:);

        dphidp(p,m) = -2 * sum(Wm(:) .* dJm(:));
    end
end

end

function dphidp = compute_cost_function_gradient_a_opt_optimized(...
    img,sensor_locations,Gamma_prior_inv,Gamma_noise_inv,A,dim)

n_sensors = numel(img.fwd_model.sensors);
n_stim = numel(img.fwd_model.stimulation);
n_elem = size(img.fwd_model.elems,1);

img = assign_sensor_locations(img,sensor_locations);

J = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,dim);
dJds = compute_dJ_xyz(img,dim);

H = J.'*Gamma_noise_inv*J + Gamma_prior_inv;

L = chol(H,'lower');

% Compute H^{-2} through Cholesky factorization and linear system solves

% Solve L'Y = I -> LH^{-1} = Y -> L'W = H^{-1} -> L H^{-2} = W
Y = L'\eye(n_elem);
X = L\Y;
W = L'\X;
inv_H2 = L\W;

W = Gamma_noise_inv*J*inv_H2;

dphidp = zeros(3,n_sensors);

for p = 1:3
    for m = 1:n_sensors
        ids = m : n_sensors : n_sensors*n_stim;

        dJm = reshape(dJds{p}(m,:,:), [n_stim, n_elem]);
        Wm  = W(ids,:);

        dphidp(p,m) = -2 * sum(Wm(:) .* dJm(:));
    end
end

end


% For 3-axis MDEIT
function dphidp = compute_cost_function_gradient_a_opt_optimized_3_axis(...
    img,sensor_locations,Gamma_prior_inv,Gamma_noise_inv,A)

n_sensors = numel(img.fwd_model.sensors);
n_stim = numel(img.fwd_model.stimulation);
n_elem = size(img.fwd_model.elems,1);

img = assign_sensor_locations(img,sensor_locations);

Jx = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,1);
Jy = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,2);
Jz = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,3);

% dJxds2 = compute_dJ_xyz(img,1);
dJ = compute_dJxyz_xyz(img);

dJxds = cell(1,3);
dJyds = cell(1,3);
dJzds = cell(1,3);

for p = 1:3
    dJxds{p} = dJ{1,p};
    dJyds{p} = dJ{2,p};
    dJzds{p} = dJ{3,p};
end

J = [Jx;Jy;Jz];

H = J.'*Gamma_noise_inv*J + Gamma_prior_inv;

L = chol(H,'lower');

% Compute H^{-2} through Cholesky factorization and linear system solves

% Solve L'Y = I -> LH^{-1} = Y -> L'W = H^{-1} -> L H^{-2} = W
Y = L'\eye(n_elem);
X = L\Y;
W = L'\X;
inv_H2 = L\W;

% THis is wrong I think
% Y = L \ (L' \ J');      % H^{-1} J'
% Z = L \ (L' \ Y);       % H^{-2} J'
% W = Gamma_noise_inv * J * Z';


W = Gamma_noise_inv*J*inv_H2;

dphidp = zeros(3,n_sensors);

block_size = n_sensors * n_stim;

for p = 1:3
    for m = 1:n_sensors
        % indices inside one block (x, y, or z)
        ids_local = m : n_sensors : block_size;

        % --- X component ---
        ids_x = ids_local;
        dJx_m = reshape(dJxds{p}(m,:,:), [n_stim, n_elem]);
        Wx_m  = W(ids_x,:);
        term_x = sum(Wx_m(:) .* dJx_m(:));

        % --- Y component ---
        ids_y = ids_local + block_size;
        dJy_m = reshape(dJyds{p}(m,:,:), [n_stim, n_elem]);
        Wy_m  = W(ids_y,:);
        term_y = sum(Wy_m(:) .* dJy_m(:));

        % --- Z component ---
        ids_z = ids_local + 2*block_size;
        dJz_m = reshape(dJzds{p}(m,:,:), [n_stim, n_elem]);
        Wz_m  = W(ids_z,:);
        term_z = sum(Wz_m(:) .* dJz_m(:));

        dphidp(p,m) = -2 * (term_x + term_y + term_z);
    end
end
end


function dphidp = compute_cost_function_gradient_d_opt_optimized_3_axis(...
    img,sensor_locations,Gamma_prior_inv,Gamma_noise_inv,A)

n_sensors = numel(img.fwd_model.sensors);
n_stim = numel(img.fwd_model.stimulation);
n_elem = size(img.fwd_model.elems,1);

img = assign_sensor_locations(img,sensor_locations);

Jx = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,1);
Jy = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,2);
Jz = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,3);

% dJxds2 = compute_dJ_xyz(img,1);
dJ = compute_dJxyz_xyz(img);

dJxds = cell(1,3);
dJyds = cell(1,3);
dJzds = cell(1,3);

for p = 1:3
    dJxds{p} = dJ{1,p};
    dJyds{p} = dJ{2,p};
    dJzds{p} = dJ{3,p};
end

J = [Jx;Jy;Jz];

H = J.'*Gamma_noise_inv*J + Gamma_prior_inv;

L = chol(H,'lower');

% Compute H^{-2} through Cholesky factorization and linear system solves

% Solve L'Y = I -> LH^{-1} = Y -> L'W = H^{-1} -> L H^{-2} = W
Y = L'\eye(n_elem);
invH = L\Y;

W = Gamma_noise_inv*J*invH;

dphidp = zeros(3,n_sensors);

block_size = n_sensors * n_stim;

for p = 1:3
    for m = 1:n_sensors
        % indices inside one block (x, y, or z)
        ids_local = m : n_sensors : block_size;

        % --- X component ---
        ids_x = ids_local;
        dJx_m = reshape(dJxds{p}(m,:,:), [n_stim, n_elem]);
        Wx_m  = W(ids_x,:);
        term_x = sum(Wx_m(:) .* dJx_m(:));

        % --- Y component ---
        ids_y = ids_local + block_size;
        dJy_m = reshape(dJyds{p}(m,:,:), [n_stim, n_elem]);
        Wy_m  = W(ids_y,:);
        term_y = sum(Wy_m(:) .* dJy_m(:));

        % --- Z component ---
        ids_z = ids_local + 2*block_size;
        dJz_m = reshape(dJzds{p}(m,:,:), [n_stim, n_elem]);
        Wz_m  = W(ids_z,:);
        term_z = sum(Wz_m(:) .* dJz_m(:));

        dphidp(p,m) = -2 * (term_x + term_y + term_z);
    end
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

