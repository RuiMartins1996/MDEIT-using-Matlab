function  img = place_anomaly(img, varargin)
if nargin == 1 %Default use
    anomaly_conductivity = 3.0;
    anomaly_radius = min(max(img.fwd_model.nodes,1)-min(img.fwd_model.nodes,1))/5;
    anomaly_center = mean(img.fwd_model.nodes,1);

    img = placeSphericalAnomaly(img,anomaly_conductivity,anomaly_radius,anomaly_center);
else %If struct modelParameters is given
    model_parameters = varargin{1};

    assert(isfield(model_parameters,'anomaly'),'No field "anomaly" found for second argument')
    assert(isfield(model_parameters.anomaly,'radius'),'No field "radius" found in anomaly');
    assert(isfield(model_parameters.anomaly,'conductivity'),'No field "conductivity" found in anomaly');

    switch model_parameters.anomaly.type
        case 'spherical'
            anomaly_radius = model_parameters.anomaly.radius;
            anomaly_conductivity = model_parameters.anomaly.conductivity;

            if isfield(model_parameters.anomaly,'position')
                anomaly_center = model_parameters.anomaly.position;

                % Check if spherical anomaly is contained on the domain.
                % ASSUME cylindrical domain!
                cylinder_center = mean(img.fwd_model.nodes,1);
                cylinder_center = cylinder_center(:);
                
                distance_sphere_cylinder_center = norm(cylinder_center(:)-anomaly_center(:),2);
                

                assert(...
                    distance_sphere_cylinder_center+anomaly_radius < model_parameters.radius &&...
                    abs(anomaly_center(3)-cylinder_center(3))+anomaly_radius<model_parameters.height/2 ...
                    ,'Anomaly will not fit at this position');
            else
                anomaly_center = findSphericalAnomalyCenter(model_parameters);
            end

            img = placeSphericalAnomaly(img,anomaly_conductivity,anomaly_radius,anomaly_center);
        case 'cylindrical'
            
            anomaly_radius = model_parameters.anomaly.radius;
            anomaly_conductivity = model_parameters.anomaly.conductivity;
            
            if isfield(model_parameters.anomaly,'position')
                anomaly_center = model_parameters.anomaly.position;
            else
                anomaly_center = mean(img.fwd_model.nodes);
            end

            img = place_cylindrical_anomaly(img,anomaly_conductivity,anomaly_radius,anomaly_center);
        otherwise
            error(strcat('Type',string(model_parameters.anomaly.type),' is not implemented'));
    end
end
end



%% FUNCTION: placeSphericalAnomaly
function img = placeSphericalAnomaly(img,anomaly_conductivity,anomaly_radius,anomaly_center)

str = strcat('(x-',num2str(anomaly_center(1)),...
    ').^2 + (y-',num2str(anomaly_center(2)),...
    ').^2 + (z-',num2str(anomaly_center(3)),...
    ').^2 <',num2str(anomaly_radius),'^2');

select_fcn = inline(str,'x','y','z');

memb_frac = elem_select(img.fwd_model, select_fcn);
img.elem_data(memb_frac~=0) = anomaly_conductivity;%*memb_frac(memb_frac~=0);

end

%% FUNCTION:
function img = place_cylindrical_anomaly(img,anomaly_conductivity,anomaly_radius,anomaly_center)

str = strcat('(x-',num2str(anomaly_center(1)),...
    ').^2 + (y-',num2str(anomaly_center(2)),...
    ').^2 <',num2str(anomaly_radius),'^2');

select_fcn = inline(str,'x','y','z');

memb_frac = elem_select(img.fwd_model, select_fcn);
img.elem_data(memb_frac~=0) = anomaly_conductivity;%*memb_frac(memb_frac~=0);

end

%% FUNCTION: findSphericalAnomalyCenter
function anomalyCenter = findSphericalAnomalyCenter(modelParameters)
% Find an admisseable center to place an anomaly of anomalyRadius

anomalyRadius = modelParameters.anomaly.radius;

numOfRings = modelParameters.numOfRings;

dh = modelParameters.height/(numOfRings+1);
ring_vert_pos = dh:dh:modelParameters.height-dh;

%Anomaly must be inside imaging plane, so ring_vert_pos must be
%accounted!
if numel(ring_vert_pos)>1
    upperBound = max(ring_vert_pos);
    lowerBound = min(ring_vert_pos);

    %Ensure that it will fit inside the imaging plane
    ub = upperBound-anomalyRadius;
    lb = lowerBound+anomalyRadius;

    assert(ub>lb,'Anomaly will not fit in this model')

else
    assert(modelParameters.height>anomalyRadius,'Anomaly will not fit the model');

    upperBound = modelParameters.height;
    lowerBound = 0;
end

%Generate random cylindrical coordinates
absAnomalyZ = lb + (ub - lb) * rand(1);

ub = modelParameters.radius-anomalyRadius;
lb = 0;

if ub<lb
    error('Anomaly will not fit in this model')
end
absAnomalyR = lb + (ub - lb) * rand(1);

absAnomalyTheta = 2*pi*rand(1);

% Convert to cartesian coordinates
absAnomalyX = absAnomalyR*cos(absAnomalyTheta);
absAnomalyY = absAnomalyR*sin(absAnomalyTheta);

absAnomalyCenter = [absAnomalyX;absAnomalyY;absAnomalyZ];

% Compute anomaly center
anomalyCenter = absAnomalyCenter;
end
