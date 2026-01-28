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
