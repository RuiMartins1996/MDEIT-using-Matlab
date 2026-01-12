classdef testRMatrices < matlab.unittest.TestCase

    properties
        testParameters;
    end

    methods(TestClassSetup)
        % Shared setup for the entire test class
    end

    methods(TestMethodSetup)
        function setup(testCase)

            % Prepare workspace
            cd("C:\Users\RuiMartins\Desktop\SVD Comparison");
            addpath(genpath("functions"));
            addpath(genpath("libraries"));

            % Setup EIDORS
            eidorsFolder = setupEidors(cd);
            clc;

            run("globalParameters.m")

            rng(1);

            % Define test parameters

            % The sensor position
            testCase.testParameters.rm = [0 0 5];

            % Current density vector
            testCase.testParameters.J = [1 1 1];

            % Error threshold
            testCase.testParameters.errorThresh = 1e-1;
        end
    end

    methods(Test)
        function testExactExpressionForRMatrices(testCase)

            %Extract Parameters
            errorThresh = testCase.testParameters.errorThresh;

            height = 1;
            radius = 1;
            maxsz = 0.5;
            numOfElectrodesPerRing = 3;
            numOfRings = 2;

            numOfSensors = 10;
            R = 5;

            % Assemble forward model
            dh = height/(numOfRings+1);
            ring_vert_pos = dh:dh:height-dh;

            %Create  forward model
            fmdl= ng_mk_cyl_models(...
                [height,radius,maxsz],[numOfElectrodesPerRing,ring_vert_pos],[0.2,0.2,maxsz]);

            options{1} = numOfSensors;
            options{2} = R;
            options{3} = sum(fmdl.nodes,1)/length(fmdl.nodes);
            sensorLocations = placeMagnetometers(options,'spherical');
            clc;

            fprintf("Test Parameters: \n"+ ...
                "Error threshold: %g (%%)\n",testCase.testParameters.errorThresh);

            %test
            [Rx,Ry,Rz] = testRMatrices.computeRmatricesNewtonCotes1(fmdl,sensorLocations);


            % Compute R matrices with two methods
            [Gx,Gy,Gz,V,elementCentroids] = computeGradientMatrix(fmdl);
            [RxCentroid,RyCentroid,RzCentroid] = computeRmatrices(fmdl,elementCentroids,sensorLocations);

            [RxAnalytic,RyAnalytic,RzAnalytic] = computeRmatricesGeneralized(fmdl,sensorLocations);

            % The result of this one should be the same as the centroid quadrature
            [RxKeast0,RyKeast0,RzKeast0] = testRMatrices.computeRmatricesKeast0(fmdl,sensorLocations);

            [RxNewtonCotes4,RyNewtonCotes4,RzNewtonCotes4] = testRMatrices.computeRmatricesNewtonCotes4(fmdl,sensorLocations);

            [RxNewtonCotes2,RyNewtonCotes2,RzNewtonCotes2] = testRMatrices.computeRmatricesNewtonCotes2(fmdl,sensorLocations);

            % DONT FORGET TO CORRECT FOR VOLUME, BECAUSE THE FIRST FUNCTION IS STILL
            % INCLUDING VOLUME IN THE DEFINITION OF R
            errors = zeros(1,length(V));
            for k = 1:length(V)
                RxCentroid(:,k) = RxCentroid(:,k)*V(k);
                RyCentroid(:,k) = RyCentroid(:,k)*V(k);
                RzCentroid(:,k) = RzCentroid(:,k)*V(k);

                errors(k) = norm(RxAnalytic(:,k)-RxCentroid(:,k));
            end

            errorsKeast0Rx = abs(RxAnalytic-RxKeast0)./abs(RxAnalytic)*100;
            errorsKeast0Rx = errorsKeast0Rx(:);

            errorsKeast0Ry = abs(RyAnalytic-RyKeast0)./abs(RyAnalytic)*100;
            errorsKeast0Ry = errorsKeast0Ry(:);

            errorsKeast0Rz = abs(RzAnalytic-RzKeast0)./abs(RzAnalytic)*100;
            errorsKeast0Rz = errorsKeast0Rz(:);


            errorsNewtonCotes2Rx = abs(RxAnalytic-RxNewtonCotes2)./abs(RxAnalytic)*100;
            errorsNewtonCotes2Rx = errorsNewtonCotes2Rx(:);

            errorsNewtonCotes2Ry = abs(RyAnalytic-RyNewtonCotes2)./abs(RyAnalytic)*100;
            errorsNewtonCotes2Ry = errorsNewtonCotes2Ry(:);

            errorsNewtonCotes2Rz = abs(RzAnalytic-RzNewtonCotes2)./abs(RzAnalytic)*100;
            errorsNewtonCotes2Rz = errorsNewtonCotes2Rz(:);


            errorsNewtonCotes4Rx = abs(RxAnalytic-RxNewtonCotes4)./abs(RxAnalytic)*100;
            errorsNewtonCotes4Rx = errorsNewtonCotes4Rx(:);

            errorsNewtonCotes4Ry = abs(RyAnalytic-RyNewtonCotes4)./abs(RyAnalytic)*100;
            errorsNewtonCotes4Ry = errorsNewtonCotes4Ry(:);

            errorsNewtonCotes4Rz = abs(RzAnalytic-RzNewtonCotes4)./abs(RzAnalytic)*100;
            errorsNewtonCotes4Rz = errorsNewtonCotes4Rz(:);

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
            fprintf('%-15s%-15.9g%-25.9g%-10.9g\n','Keast0',mean(errorsKeast0Rx),std(errorsKeast0Rx),max(errorsKeast0Rx));
            fprintf('%-15s%-15.9g%-25.9g%-10.9g\n','NewtonCotes4',mean(errorsNewtonCotes2Rx),std(errorsNewtonCotes2Rx),max(errorsNewtonCotes2Rx));
            fprintf('%-15s%-15.9g%-25.9g%-10.9g\n','NewtonCotes35',mean(errorsNewtonCotes4Rx),std(errorsNewtonCotes4Rx),max(errorsNewtonCotes4Rx));

            fprintf('%s\n', repmat('-', 1, 30));

            fprintf('%-15s\n','Ry:');
            fprintf('%-15s%-15s%-25s%-10s\n','Quadrature', 'Mean (%)', 'Standard Deviation (%)','Maximum (%)');
            fprintf('%-15s%-15.9g%-25.9g%-10.9g\n','Keast0',mean(errorsKeast0Ry),std(errorsKeast0Ry),max(errorsKeast0Ry));
            fprintf('%-15s%-15.9g%-25.9g%-10.9g\n','NewtonCotes4',mean(errorsNewtonCotes2Ry),std(errorsNewtonCotes2Ry),max(errorsNewtonCotes2Ry));
            fprintf('%-15s%-15.9g%-25.9g%-10.9g\n','NewtonCotes35',mean(errorsNewtonCotes4Ry),std(errorsNewtonCotes4Ry),max(errorsNewtonCotes4Ry));
            fprintf('%s\n', repmat('-', 1, 30));

            fprintf('%-15s\n','Rz:');
            fprintf('%-15s%-15s%-25s%-10s\n','Quadrature', 'Mean (%)', 'Standard Deviation (%)','Maximum (%)');
            fprintf('%-15s%-15.9g%-25.9g%-10.9g\n','Keast0',mean(errorsKeast0Rz),std(errorsKeast0Rz),max(errorsKeast0Rz));
            fprintf('%-15s%-15.9g%-25.9g%-10.9g\n','NewtonCotes4',mean(errorsNewtonCotes2Rz),std(errorsNewtonCotes2Rz),max(errorsNewtonCotes2Rz));
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
