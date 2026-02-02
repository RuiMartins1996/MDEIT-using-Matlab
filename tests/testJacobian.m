classdef testJacobian < matlab.unittest.TestCase
    
    methods(TestClassSetup)
        % Shared setup for the entire test class
    end
    
    methods(TestMethodSetup)
        % Setup for each test
    end
    
    methods(Test)
        % Test methods
        
        % Test the jacobian computation when the measurement axis is the
        % default cartesian
        function testJacobianComputationCartesian(testCase)
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

            current_amplitude = 2.4e-3/I0;
            background_conductivity = 3.28e-1/sigma0;

            
            model_parameters = create_default_3d_model_parameters(l0, z0, sigma0, I0);
            
            model_parameters.maxsz = min(model_parameters.height,model_parameters.radius);
            model_parameters.height =1;
            model_parameters.numOfRings = 1;
            model_parameters.numOfElectrodesPerRing = 4;
            model_parameters.numOfSensors = 2;
            model_parameters.sensorRadius = model_parameters.radius*2;
            model_parameters.material = struct();

            %% Create a very coarse forward model
            [~,fmdls] = mk_mdeit_model(model_parameters,model_folder);
            fmdl = fmdls{1};
            %% Stimulation pattern is not the default. For now, edit it manually
            inj = [0 3]; %skip 2 pattern (pg 172)
            meas = [0 3]; %for EIT, skip2 measurement protocol was used
            stimulation = mk_stim_patterns(model_parameters.numOfElectrodesPerRing,model_parameters.numOfRings,inj,meas,{},current_amplitude);

            fmdl.stimulation = stimulation;

            %Make homogeneous image
            imgh = mk_image_mdeit(fmdl,background_conductivity);
            
            %% Compute jacobian
            verbose = true;
            lambdatimesdAdp = []; %legacy
            
            A = @(x) M(imgh,x);
            recon_mode = 'mdeit1';

            select_sensor_axis = 1;
            J_mdeit_1 = calc_jacobian_mdeit(imgh,imgh.elem_data,lambdatimesdAdp,A,recon_mode,select_sensor_axis,verbose);

            select_sensor_axis = 2;
            J_mdeit_2 = calc_jacobian_mdeit(imgh,imgh.elem_data,lambdatimesdAdp,A,recon_mode,select_sensor_axis,verbose);

            select_sensor_axis = 3;
            J_mdeit_3 = calc_jacobian_mdeit(imgh,imgh.elem_data,lambdatimesdAdp,A,recon_mode,select_sensor_axis,verbose);
            
            recon_mode = 'mdeit3';
            J_mdeit_3d = calc_jacobian_mdeit(imgh,imgh.elem_data,lambdatimesdAdp,A,recon_mode);

            % Compute jacobian with finite differences
            [J_mdeit_fd_1,J_mdeit_fd_2,J_mdeit_fd_3] = calc_jacobian_finite_differences(imgh,min(imgh.elem_data)*1e-6);

            error_1_percentage = abs(J_mdeit_fd_1(:)-J_mdeit_1(:))./abs(J_mdeit_1(:))*100;
            error_2_percentage = abs(J_mdeit_fd_2(:)-J_mdeit_2(:))./abs(J_mdeit_2(:))*100;
            error_3_percentage = abs(J_mdeit_fd_3(:)-J_mdeit_3(:))./abs(J_mdeit_3(:))*100;
            
            n_meas = size(J_mdeit_1,1);

            num_1 = J_mdeit_1 - J_mdeit_3d(1:n_meas,:);
            den_1 = J_mdeit_3d(1:n_meas,:);

            num_2 = J_mdeit_2 - J_mdeit_3d(n_meas+1:2*n_meas,:);
            den_2 = J_mdeit_3d(n_meas+1:2*n_meas,:);

            num_3 = J_mdeit_3 - J_mdeit_3d(2*n_meas+1:3*n_meas,:);
            den_3 = J_mdeit_3d(2*n_meas+1:3*n_meas,:);

            % num_1 = J_mdeit_fd_1 - J_mdeit_3d(1:n_meas,:);
            % den_1 = J_mdeit_3d(1:n_meas,:);
            %
            % num_2 = J_mdeit_fd_2 - J_mdeit_3d(n_meas+1:2*n_meas,:);
            % den_2 = J_mdeit_3d(n_meas+1:2*n_meas,:);
            % 
            % num_3 = J_mdeit_fd_3 - J_mdeit_3d(2*n_meas+1:3*n_meas,:);
            % den_3 = J_mdeit_3d(2*n_meas+1:3*n_meas,:);

            error_3d_percentage_1 = abs(num_1(:))./abs(den_1(:))*100;
            error_3d_percentage_2 = abs(num_2(:))./abs(den_2(:))*100;
            error_3d_percentage_3 = abs(num_3(:))./abs(den_3(:))*100;
            
            testCase.verifyTrue(norm(error_1_percentage,'inf')<1);
            testCase.verifyTrue(norm(error_2_percentage,'inf')<1);
            testCase.verifyTrue(norm(error_3_percentage,'inf')<1);

            testCase.verifyTrue(norm(error_3d_percentage_1,'inf')<1);
            testCase.verifyTrue(norm(error_3d_percentage_2,'inf')<1);
            testCase.verifyTrue(norm(error_3d_percentage_3,'inf')<1);

        end
        
        % Test the jacobian computation when the measuremnt axis is
        % cylindrical
        function testJacobianComputationCylindrical(testCase)
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

            current_amplitude = 2.4e-3/I0;
            background_conductivity = 3.28e-1/sigma0;

            model_parameters = create_default_3d_model_parameters(l0, z0, sigma0, I0);
            
            model_parameters.maxsz = min(model_parameters.height,model_parameters.radius);
            model_parameters.height =1;
            model_parameters.numOfRings = 1;
            model_parameters.numOfElectrodesPerRing = 4;
            model_parameters.numOfSensors = 4;
            model_parameters.sensorRadius = model_parameters.radius*2;
            model_parameters.material = struct();
            model_parameters.measurementAxisType = 'cylindrical';

            %% Create a very coarse forward model
            [model_parameters,fmdls] = mk_mdeit_model(model_parameters,model_folder);
            fmdl = fmdls{1};
            %% Stimulation pattern is not the default. For now, edit it manually
            inj = [0 3]; %skip 2 pattern (pg 172)
            meas = [0 3]; %for EIT, skip2 measurement protocol was used
            stimulation = mk_stim_patterns(model_parameters.numOfElectrodesPerRing,model_parameters.numOfRings,inj,meas,{},current_amplitude);

            fmdl.stimulation = stimulation;

            %Make homogeneous image
            imgh = mk_image_mdeit(fmdl,background_conductivity);
            
            % %% 
            % hold on
            % show_fem(fmdl);
            % plot_sensors(fmdl);

            %% Compute jacobian

            verbose = true;
            lambdatimesdAdp = []; %legacy

            A = @(x) M(imgh,x);
            recon_mode = 'mdeit1';

            select_sensor_axis = 1;
            J_mdeit_1 = calc_jacobian_mdeit(imgh,imgh.elem_data,lambdatimesdAdp,A,recon_mode,select_sensor_axis,verbose);
            
            % Compute jacobian with finite-differences
            [J_mdeit_fd_1,J_mdeit_fd_2,J_mdeit_fd_3] = calc_jacobian_finite_differences(imgh,min(imgh.elem_data)*1e-6);
            
            error_1_percentage = abs(J_mdeit_fd_1(:)-J_mdeit_1(:))./abs(J_mdeit_1(:))*100;
            % error_2_percentage = abs(J_mdeit_fd_2(:)-J_mdeit_2(:))./abs(J_mdeit_2(:))*100;
            % error_3_percentage = abs(J_mdeit_fd_3(:)-J_mdeit_3(:))./abs(J_mdeit_3(:))*100;

        end
    end
    
end



function [J1,J2,J3] = calc_jacobian_finite_differences(img,delta)
    
    n_sensors = numel(img.fwd_model.sensors);
    n_inj = numel(img.fwd_model.stimulation);
    
    n_elem = size(img.fwd_model.elems,1);

    J1 = zeros(n_sensors*n_inj,n_elem);
    J2 = zeros(n_sensors*n_inj,n_elem);
    J3 = zeros(n_sensors*n_inj,n_elem);
    
    % Solve forward model once for unperturbed case
    data = fwd_solve_mdeit(img);
    B1 = data.Bx;
    B2 = data.By;
    B3 = data.Bz;
    
    elem_data = img.elem_data;

    for k = 1:n_elem
        fprintf('Element %i of %i \n',k,n_elem)
        elem_data_k = elem_data;
        elem_data_k(k) = elem_data(k)+delta;
        
        img.elem_data = elem_data_k;

        data = fwd_solve_mdeit(img);
        
        dB1_k = -(B1(:)-data.Bx(:))/delta;
        dB2_k = -(B2(:)-data.By(:))/delta;
        dB3_k = -(B3(:)-data.Bz(:))/delta;

        J1(:,k) = dB1_k;
        J2(:,k) = dB2_k;
        J3(:,k) = dB3_k;
    end
end
%% FUNCTIONS
function out = M(img,sigma)

numNodes = size(img.fwd_model.nodes,1);

img.elem_data = sigma;
s_mat = system_mat_1st_order(img);

Ac = s_mat.E(1:numNodes,1:numNodes);
Ae = s_mat.E(1:numNodes,numNodes+1:end);
Ad = s_mat.E(numNodes+1:end,numNodes+1:end);

out = Ac-Ae*inv(Ad)*Ae';
end
