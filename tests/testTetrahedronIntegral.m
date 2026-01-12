classdef testTetrahedronIntegral < matlab.unittest.TestCase

    %% (Optional) Properties shared by tests
    properties
        testParameters  % You can define any data needed across test methods
    end
    %% Code that runs before each test method
    methods (TestMethodSetup)
        function setup(testCase)
            clc;

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
            testCase.testParameters.R = 10;
            testCase.testParameters.numOfTetrahedra = 10;
            testCase.testParameters.numOfSensors = 10;
            testCase.testParameters.errorThresh = 1e-3;
            
            % Display parameters
            fprintf("Test Parameters: \n" + ...
                "R: %g\n" + ...
                "Number of tetrahedra: %i\n" + ...
                "Number of sensors: %i\n" + ...
                "Error threshold: %g (%%)\n",...
                testCase.testParameters.R,...
                testCase.testParameters.numOfTetrahedra,...
                testCase.testParameters.numOfSensors,...
                testCase.testParameters.errorThresh);
        end
    end

    %% Code that runs after each test method
    methods (TestMethodTeardown)
        function teardown(testCase)

        end
    end

    %% Code that runs once before all test methods
    methods (TestClassSetup)
        function classSetup(testCase)



        end
    end

    %% Code that runs once after all test methods
    methods (TestClassTeardown)
        function classTeardown(testCase)

        end
    end

    %% Test Methods
    methods (Test)
        %% TEST: Compares exact integration to the numerical integration on the tetrahedron
        function testIfAbsoluteErrorIsLessThanThreshold(testCase)
            %% Fetch parameters
            R = testCase.testParameters.R;
            numOfTetrahedra = testCase.testParameters.numOfTetrahedra;
            numOfSensors = testCase.testParameters.numOfSensors;
            errorThresh = testCase.testParameters.errorThresh;

            %% Create sensor locations
            v = randn(3,numOfSensors);
            normv = sqrt(sum(v.^2,1));

            rms = zeros(size(v));
            for i = 1:size(v,2)
                rms(:,i) = 2*R*v(:,i)./normv(i);
            end

            %% Generate random tetrahedra on a sphere of radius R and center (0,0,0)

            tetrahedra = zeros(3,4,numOfTetrahedra);

            for i = 1:numOfTetrahedra
                tet = zeros(3,4);
                for j = 1:4
                    theta = 2*pi*rand();
                    phi = acos(2*rand()-1);
                    tet(:,j) = [R*cos(theta)*sin(phi), ...
                        R*sin(theta)*sin(phi), ...
                        R*cos(phi)]';
                end
                tetrahedra(:,:,i) = tet;
            end
            
            %% --- Helper function ---
            function val = integrand_mapped(r, s, t, V1, T,f)
                % r, s, t can be arrays
                rst = cat(1, r(:)', s(:)', t(:)');  % size: 3 x N
                xyz = V1' + T * rst;               % 3 x N
                x = reshape(xyz(1,:), size(r));
                y = reshape(xyz(2,:), size(r));
                z = reshape(xyz(3,:), size(r));

                val = f(x,y,z); % Replace this with your own function
            end
            %% Loop over tetrahedra and compute the integral

            for t = 1:numOfTetrahedra
                vertices = tetrahedra(:,:,t)';

                V1 = vertices(1,:);
                V2 = vertices(2,:);
                V3 = vertices(3,:);
                V4 = vertices(4,:);
                % Create EIT-FEM model
                fmdl = struct();
                fmdl.type = 'fwd_model';
                fmdl.nodes(1:4,:) = vertices;
                fmdl.elems(1,1:4) = [1 2 3 4];

                fmdl = fix_model(fmdl);

                %% Generate random unitary current density vector
                
                a = -1;
                b = 1;
                J = a + (b-a).*rand(1,3);
                J = J./norm(J);
                
                %% Compare the analytic vs numerical integral on a single tetrahedron and single sensor
                
                for j = 1:numOfSensors
                    
                    rm = rms(:,j)';

                    % Exact integral
                    Ikm = zeros(1,3);
                    [L,n] = computeIntegralAndNormalsInTetrahedron(vertices,rm);

                    for i = 1:4
                        Ikm = Ikm + cross(n(:,i),J)*L(i);
                    end

                    Ikm = -Ikm;

                    % Numerical integral using matlab builtin integral3

                    T = [V2 - V1; V3 - V1; V4 - V1]'; % 3x3 matrix
                    detJ = abs(det(T));               % Volume scale factor

                    % Define function f(x,y,z) to integrate
                    fx = @(x,y,z) (J(2)*(rm(3)-z)-J(3)*(rm(2)-y)) ./ (sqrt((rm(1)-x).^2 + (rm(2)-y).^2 + (rm(3)-z).^2).^3);
                    fy = @(x,y,z) -(J(1)*(rm(3)-z)-J(3)*(rm(1)-x)) ./ (sqrt((rm(1)-x).^2 + (rm(2)-y).^2 + (rm(3)-z).^2).^3);
                    fz = @(x,y,z) (J(1)*(rm(2)-y)-J(2)*(rm(1)-x)) ./ (sqrt((rm(1)-x).^2 + (rm(2)-y).^2 + (rm(3)-z).^2).^3);

                    f = {fx,fy,fz};

                    I = zeros(1,3);
                    for i = 1:3
                        % Define the integrand as a function of r, s, t (reference tetrahedron)
                        f_rst = @(r,s,t) ...
                            integrand_mapped(r, s, t, V1, T,f{i});

                        rmin = 0;
                        rmax = 1;
                        smin = @(r) zeros(size(r));      % vectorized zeros
                        smax = @(r) 1 - r;
                        tmin = @(r,s) zeros(size(r));    % vectorized zeros
                        tmax = @(r,s) 1 - r - s;

                        % Perform the integral over reference tetrahedron (r,s,t)
                        I_ref = integral3( ...
                            @(r,s,t) f_rst(r,s,t), ...
                            rmin,rmax, ...
                            smin, smax, ...
                            tmin, tmax);

                        % Scale by |detJ| (Jacobian determinant)
                        I(i) = detJ * I_ref;
                    end
                    
                    %% Test
                    testCase.verifyTrue(max(abs(I-Ikm))<errorThresh,...
                        sprintf('The difference between analytical and numerical integral on single tetrahedron is less than %d\n',errorThresh));
                end

            end
            %% Present results for a single tetrahedron
            fprintf('%s\n', repmat('-', 1, 30));

            aError = abs((Ikm-I));

            fprintf('%-10s\n', 'Relative Error for last tetrahedron tested:');
            fprintf('%-20s %-15s %-20s %-20s\n','', 'Exact Integral', 'Matlab Quadrature','Absolute Error');
            fprintf('%-20s %-15.4g %-20.4g %-15.2d\n','x - component',Ikm(1),I(1),aError(1));
            fprintf('%-20s %-15.4g %-20.4g %-15.2d\n','y - component',Ikm(2),I(2),aError(2));
            fprintf('%-20s %-15.4g %-20.4g %-15.2d\n','z - component',Ikm(3),I(3),aError(3));

            fprintf('%s\n', repmat('-', 1, 30));
        end
    end

    %% Static methods
    methods (Static, Access = private)

    end
end
