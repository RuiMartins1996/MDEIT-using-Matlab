classdef testTriangleIntegral < matlab.unittest.TestCase

    %% (Optional) Properties shared by tests
    properties
        testParameters  % You can define any data needed across test methods
    end

    %% (Optional) Setup and Teardown Methods
    methods (TestMethodSetup) % Code that runs before each test method
        function setup(testCase)
            clc;

            % Prepare workspace
            cd("C:\Users\RuiMartins\Desktop\SVD Comparison");
            addpath(genpath("functions"));

            rng(1);
            
            % Define test parameters
            testCase.testParameters.R = 10;
            testCase.testParameters.numOfTriangles = 1000;
            testCase.testParameters.Nsensors = 100;
            testCase.testParameters.errorThresh = 1;

            %% Display the testCase parameters
            R = testCase.testParameters.R;
            numOfTriangles = testCase.testParameters.numOfTriangles;
            Nsensors = testCase.testParameters.Nsensors;

            errorThresh = testCase.testParameters.errorThresh;

            fprintf("Test Parameters: \n" + ...
                "R: %g\n" + ...
                "Number of triangles: %i\n" + ...
                "Number of sensors: %i\n" + ...
                "Error threshold: %g (%%)\n",R,numOfTriangles,Nsensors,errorThresh);
            fprintf('%s\n', repmat('-', 1, 30));
        end
    end

    methods (TestMethodTeardown) % Code that runs after each test method
        function teardown(testCase)

        end
    end

    methods (TestClassSetup) % Code that runs once before all test methods
        function classSetup(testCase)

        end
    end

    methods (TestClassTeardown) % Code that runs once after all test methods
        function classTeardown(testCase)

        end
    end

    %% Test Methods
    methods (Test)
        %% TEST: Tests if quadrature integration of the triangle integral is close to the analytical expression
        function testIfRelativeErrorsAreLessThanThreshold(testCase) % A test case that uses the shared property
            %% Fetch test parameters and display
            
            R = testCase.testParameters.R;
            numOfTriangles = testCase.testParameters.numOfTriangles;
            Nsensors = testCase.testParameters.Nsensors;

            errorThresh = testCase.testParameters.errorThresh;
            
            %% Create sensor locations
            v = randn(3,Nsensors);
            normv = sqrt(sum(v.^2,1));

            rm = zeros(size(v));
            for i = 1:size(v,2)
                rm(:,i) = 2*R*v(:,i)./normv(i);
            end

            %% Generate random triangles on a sphere of radius R and center (0,0,0)

            triangles = zeros(3,3,numOfTriangles);

            for i = 1:numOfTriangles/2

                tri = zeros(3,3);

                for j = 1:3
                    theta = 2*pi*rand();
                    phi = acos(2*rand()-1);

                    tri(:,j) = [R*cos(theta)*sin(phi),R*sin(theta)*sin(phi),R*cos(phi)]';
                end

                triangles(:,:,i) = tri;
            end

            %% Generate random right triangles

            for i = numOfTriangles/2+1:numOfTriangles
                tri = zeros(3,3);

                theta = 2*pi*rand();
                phi = acos(2*rand()-1);

                tri(:,1) = [R*cos(theta)*sin(phi),R*sin(theta)*sin(phi),R*cos(phi)]';

                [A,B] = testTriangleIntegral.find_diameter_points_for_right_triangle(tri(:,1), R);

                tri(:,2) = A;
                tri(:,3) = B;
                
                triangles(:,:,i) = tri;
            end

            %% Integrate with multiple quadrature rules

            integralAnalytical = zeros(size(rm,2),numOfTriangles);
            integralNumericalDunavant = zeros(size(rm,2),numOfTriangles);
            integralNumerical7pts =  zeros(size(rm,2),numOfTriangles);
            integralNumericalCentroid = zeros(size(rm,2),numOfTriangles);
             
            for n = 1:size(rm,2)
                for i = 1:numOfTriangles

                    r = rm(:,n);

                    v1 = triangles(:,1,i);
                    v2 = triangles(:,2,i);
                    v3 = triangles(:,3,i);

                    integralAnalytical(n,i) = computeIntegralInTriangle(r,v1,v2,v3);
                    
                    integralNumericalDunavant(n,i) = testTriangleIntegral.integrateOverTriangleDunavant(v1, v2, v3, r);

                    integralNumerical7pts(n,i) = testTriangleIntegral.integrateOverTriangle7pts(v1, v2, v3, r);

                    integralNumericalCentroid(n,i) = testTriangleIntegral.integrateOverTriangleCentroid(v1, v2, v3, r);
                end
            end
            
            %% Compute relative errors
            rError7pts = abs(integralNumerical7pts(:)-integralAnalytical(:))./abs(integralAnalytical(:))*100;
            rErrorCentroid = abs(integralNumericalCentroid(:)-integralAnalytical(:))./abs(integralAnalytical(:))*100;
            rErrorDunavant = abs(integralNumericalDunavant(:)-integralAnalytical(:))./abs(integralAnalytical(:))*100;
            
            %% Test
            testCase.verifyTrue(all(rErrorCentroid < 10*errorThresh),sprintf('Not all centroid quadrature relative errors are less than %g (%%)', 10*errorThresh));
            testCase.assertTrue(all(rError7pts < errorThresh),sprintf('Not all 7pts quadrature relative errors are less than %g (%%)', errorThresh));
            testCase.assertTrue(all(rErrorDunavant < errorThresh),sprintf('Not all Dunavant quadrature relative errors are less than %g (%%)', errorThresh));
            
            %% Present results
            fprintf('%-10s\n', 'Relative Error:');
            fprintf('%-10s %-15s %-20s %-15s\n', 'Quadrature', 'Average (%)', 'Standard Deviation (%)', 'Maximum (%)');
            fprintf('%-10s %-15.4g %-20.4g %-15.4g\n','Centroid', ...
                mean(rErrorCentroid), std(rErrorCentroid), max(rErrorCentroid));
            fprintf('%-10s %-15.4g %-20.4g %-15.4g\n','7pts', ...
                mean(rError7pts), std(rError7pts), max(rError7pts));
            fprintf('%-10s %-15.4g %-20.4g %-15.4g\n','Dunavant', ...
                mean(rErrorDunavant), std(rErrorDunavant), max(rErrorDunavant));
            
            fprintf('%s\n', repmat('-', 1, 30));

        end
    end
    
    %% Static methods
    methods (Static, Access = private)

        function I = integrateOverTriangleDunavant(v1, v2, v3, rm)
            % Ensure all vectors are column vectors
            v1 = v1(:); v2 = v2(:); v3 = v3(:); rm = rm(:);

            % Triangle area
            normal = cross(v2 - v1, v3 - v1);
            area = 0.5 * norm(normal);

            bary=[
                3.33333333333333314829616256247391e-01 3.33333333333333314829616256247391e-01 3.33333333333333314829616256247391e-01 9.71357962827994192434033493555035e-02
                4.89682519198738008814331124085584e-01 4.89682519198738008814331124085584e-01 2.06349616025239823713377518288326e-02 3.13347002271391339434103429084644e-02
                4.89682519198738008814331124085584e-01 2.06349616025239823713377518288326e-02 4.89682519198738008814331124085584e-01 3.13347002271391339434103429084644e-02
                2.06349616025239823713377518288326e-02 4.89682519198738008814331124085584e-01 4.89682519198738008814331124085584e-01 3.13347002271391339434103429084644e-02
                4.37089591492937024064246998023009e-01 4.37089591492937024064246998023009e-01 1.25820817014125951871506003953982e-01 7.78275410047743337882408809491608e-02
                4.37089591492937024064246998023009e-01 1.25820817014125951871506003953982e-01 4.37089591492937024064246998023009e-01 7.78275410047743337882408809491608e-02
                1.25820817014125951871506003953982e-01 4.37089591492937024064246998023009e-01 4.37089591492937024064246998023009e-01 7.78275410047743337882408809491608e-02
                1.88203535619032996661914580727171e-01 1.88203535619032996661914580727171e-01 6.23592928761933951165019607287832e-01 7.96477389272103319939333232468925e-02
                1.88203535619032996661914580727171e-01 6.23592928761933951165019607287832e-01 1.88203535619032996661914580727171e-01 7.96477389272103319939333232468925e-02
                6.23592928761933951165019607287832e-01 1.88203535619032996661914580727171e-01 1.88203535619032996661914580727171e-01 7.96477389272103319939333232468925e-02
                4.47295133944529965663861048597028e-02 4.47295133944529965663861048597028e-02 9.10540973211094062378379021538422e-01 2.55776756586981075802800233987000e-02
                4.47295133944529965663861048597028e-02 9.10540973211094062378379021538422e-01 4.47295133944529965663861048597028e-02 2.55776756586981075802800233987000e-02
                9.10540973211094062378379021538422e-01 4.47295133944529965663861048597028e-02 4.47295133944529965663861048597028e-02 2.55776756586981075802800233987000e-02
                3.68384120547360013886439844554843e-02 2.21962989160766011043079970477265e-01 7.41198598784498008384957756788936e-01 4.32835393772891818819914533378324e-02
                3.68384120547360013886439844554843e-02 7.41198598784498008384957756788936e-01 2.21962989160766011043079970477265e-01 4.32835393772891818819914533378324e-02
                2.21962989160766011043079970477265e-01 3.68384120547360013886439844554843e-02 7.41198598784498008384957756788936e-01 4.32835393772891818819914533378324e-02
                2.21962989160766011043079970477265e-01 7.41198598784498008384957756788936e-01 3.68384120547360013886439844554843e-02 4.32835393772891818819914533378324e-02
                7.41198598784498008384957756788936e-01 3.68384120547360013886439844554843e-02 2.21962989160766011043079970477265e-01 4.32835393772891818819914533378324e-02
                7.41198598784498008384957756788936e-01 2.21962989160766011043079970477265e-01 3.68384120547360013886439844554843e-02 4.32835393772891818819914533378324e-02
                ];


            I = 0;
            for i = 1:size(bary, 1)
                l1 = bary(i, 1); l2 = bary(i, 2); l3 = bary(i, 3);
                w = bary(i, 4);

                % Map barycentric coords to physical triangle
                r = l1 * v1 + l2 * v2 + l3 * v3;

                % Evaluate integrand
                integrand = 1 / norm(r - rm);

                % Accumulate integral (area weighting comes after the loop)
                I = I + w * integrand;
            end

            % Scale by the area of the triangle
            I = I * area;
        end

        function I = integrateOverTriangle7pts(v1, v2, v3, rm)
            % Ensure all vectors are column vectors
            v1 = v1(:); v2 = v2(:); v3 = v3(:); rm = rm(:);

            % Triangle area
            normal = cross(v2 - v1, v3 - v1);
            area = 0.5 * norm(normal);

            bary=[
                3.33333333333333314829616256247391e-01 3.33333333333333314829616256247391e-01 3.33333333333333314829616256247391e-01 4.50000000000000011102230246251565e-01
                0.00000000000000000000000000000000e+00 0.00000000000000000000000000000000e+00 1.00000000000000000000000000000000e+00 5.00000000000000027755575615628914e-02
                0.00000000000000000000000000000000e+00 1.00000000000000000000000000000000e+00 0.00000000000000000000000000000000e+00 5.00000000000000027755575615628914e-02
                1.00000000000000000000000000000000e+00 0.00000000000000000000000000000000e+00 0.00000000000000000000000000000000e+00 5.00000000000000027755575615628914e-02
                5.00000000000000000000000000000000e-01 5.00000000000000000000000000000000e-01 0.00000000000000000000000000000000e+00 1.33333333333333331482961625624739e-01
                5.00000000000000000000000000000000e-01 0.00000000000000000000000000000000e+00 5.00000000000000000000000000000000e-01 1.33333333333333331482961625624739e-01
                0.00000000000000000000000000000000e+00 5.00000000000000000000000000000000e-01 5.00000000000000000000000000000000e-01 1.33333333333333331482961625624739e-01
                ];


            I = 0;
            for i = 1:size(bary, 1)
                l1 = bary(i, 1); l2 = bary(i, 2); l3 = bary(i, 3);
                w = bary(i, 4);

                % Map barycentric coords to physical triangle
                r = l1 * v1 + l2 * v2 + l3 * v3;

                % Evaluate integrand
                integrand = 1 / norm(r - rm);

                % Accumulate integral (area weighting comes after the loop)
                I = I + w * integrand;
            end

            % Scale by the area of the triangle
            I = I * area;
        end

        function I = integrateOverTriangleCentroid(v1, v2, v3, rm)
            % Ensure column vectors
            v1 = v1(:); v2 = v2(:); v3 = v3(:); rm = rm(:);

            % Compute the centroid of the triangle
            centroid = (v1 + v2 + v3) / 3;

            % Evaluate the integrand at the centroid
            integrand = 1 / norm(centroid - rm);

            % Compute the area of the triangle
            normal = cross(v2 - v1, v3 - v1);
            area = 0.5 * norm(normal);

            % Approximate the integral as area * value at centroid
            I = area * integrand;
        end

        function C = circle_center_from_sphere_plane(O, R, P, n)
            % Inputs:
            %   O - center of sphere (3x1)
            %   R - radius of sphere (scalar)
            %   P - a point on the plane (3x1)
            %   n - normal vector of plane (3x1)
            % Output:
            %   C - center of the intersection circle (3x1)

            n = n(:);  % ensure column vector
            O = O(:);
            P = P(:);

            % Vector from plane point to sphere center
            OP = O - P;

            % Project O onto the plane
            C = O - (dot(OP, n) / dot(n, n)) * n;
        end

        function [A, B] = find_diameter_points_for_right_triangle(C, R)
            % Given a point C on the surface of a sphere of radius R centered at origin,
            % find two points A and B on the sphere such that triangle ABC is right-angled at C.

            if nargin < 2
                R = norm(C);  % Infer radius from C
            end

            n = C/norm(C);

            % Get orthonormal basis (v1, v2) for the plane perpendicular to n
            vi = null(n');

            % Pick any of these two vectors to be the normal vector of a
            % perpendicular plane
            v = vi(:,randi(2));

            % Get orthonormal basis for the plane perpendicular to v
            b = null(v');

            circle_center = testTriangleIntegral.circle_center_from_sphere_plane([0,0,0], R, C, v);

            % Choose a radius vector on that plane — i.e., direction for A/B
            angle = rand() * 2 * pi;
            dir = cos(angle) * b(:,1) + sin(angle) * b(:,2);

            % A and B are endpoints of a diameter centered at C in plane orthogonal to OC
            A = circle_center + R * dir;
            B = circle_center - R * dir;
        end
    end
end
