function [sensorLocations,sensorAxes] = placeMagnetometers(options,configurationType,measurementAxisType)

if nargin < 2
    configurationType = 'cylindrical';
    measurementAxisType = 'cartesian';
end

if nargin<3
    measurementAxisType = 'cartesian';
end


switch configurationType
    case 'cylindrical'
        
        %Parse options
        if length(options)~=4 || not(iscell(options))
            error('options must be a 1x4 cell with numMagnetsPerRing,numOfRings,heights,sensorRadius');
        end
        numMagnetPerRing = options{1}; 
        numOfRings = options{2};
        heights = options{3};
        sensorRadius = options{4};
        
        %Initialize sensor locations
        numOfSensors = numMagnetPerRing*numOfRings;
        sensorLocations = zeros(numOfSensors,3);

        for m = 1:numOfRings
            for n = 1:numMagnetPerRing

                theta = n*2*pi/numMagnetPerRing;
                sensorId = n+(m-1)*numMagnetPerRing;

                sensorLocations(sensorId,1) = sensorRadius*cos(theta);
                sensorLocations(sensorId,2) = sensorRadius*sin(theta);
                sensorLocations(sensorId,3) = heights(m);
            end
        end
        
        sensorAxes = repmat(struct('axis1', [], 'axis2', [],'axis3',[]), 1, numOfSensors);
        
        for m = 1:numOfSensors
            rhat = [sensorLocations(m,1),sensorLocations(m,2),0];
            rhat = rhat/norm(rhat);
            
            thetahat = [-sensorLocations(m,2),sensorLocations(m,1),0];
            thetahat = thetahat/norm(thetahat);

            zhat = [0,0,1];

            % sensorAxes(m).axis1 = rhat;
            % sensorAxes(m).axis2 = thetahat;
            % sensorAxes(m).axis3 = zhat;

            sensorAxes(m).axis1 = [1 0 0];
            sensorAxes(m).axis2 = [0 1 0];
            sensorAxes(m).axis3 = [0 0 1];
        end

    case 'spherical'
        
        %Parse options
        if length(options)~=3 || not(iscell(options))
            error('options must be a 1x3 cell with numMagnetsPerRing,numOfRings,heights');
        end
        numOfSensors = options{1}; 
        R = options{2}; 
        center = options{3};  
        
        % Sample uniform points in the surface of the sphere
        [V,~]=SpiralSampleSphere(numOfSensors);
        
        sensorLocations = R*V+center;
        
        sensorAxes = repmat(struct('axis1', [], 'axis2', [],'axis3',[]), 1, numOfSensors);
        for m = 1:numOfSensors
            
            x = sensorLocations(m,1); y = sensorLocations(m,2); z = sensorLocations(m,3);

            r = sqrt(x^2+y^2+z^2);
            theta = acos(z/sqrt(r^2));
            phi = sign(y)*acos(x/sqrt(x^2+y^2));
            
            rhat = sin(theta)*cos(phi)*[1,0,0]+sin(theta)*sin(phi)*[0,1,0]+cos(theta)*[0,0,1];
            thetahat = cos(theta)*cos(phi)*[1,0,0]+cos(theta)*sin(phi)*[0,1,0]-sin(theta)*[0,0,1];
            phihat = -sin(phi)*[1,0,0]+cos(phi)*[0,1,0];

            % sensorAxes(m).axis1 = rhat;
            % sensorAxes(m).axis2 = thetahat;
            % sensorAxes(m).axis3 = phihat;

            sensorAxes(m).axis1 = [1 0 0];
            sensorAxes(m).axis2 = [0 1 0];
            sensorAxes(m).axis3 = [0 0 1];
        end
end








end

