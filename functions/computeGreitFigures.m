function out = computeGreitFigures(voxelImg,backgroundConductivity)

Xv = voxelImg.Xv;Yv = voxelImg.Yv;Zv = voxelImg.Zv;uVoxel = voxelImg.uVoxel;

centerOfMass = computeCenterOfMass(Xv, Yv, Zv,uVoxel,backgroundConductivity);
averagePosition = computeAveragePosition(Xv,Yv,Zv,uVoxel,backgroundConductivity);
greitResolution = computeGreitResolution(Xv,Yv,Zv,uVoxel,backgroundConductivity);
greitShapeDeformation = computeGreitShapeDeformation(Xv,Yv,Zv,uVoxel,backgroundConductivity);
greitRinging = computeGreitRinging(Xv,Yv,Zv,uVoxel,backgroundConductivity);

out = struct(...
    'centerOfMass',centerOfMass,...
    'averagePosition',averagePosition,...
    'greitResolution',greitResolution,...
    'greitShapeDeformation',greitShapeDeformation,...
    'greitRinging',greitRinging);

end


%% Functions to compute GREIT targets
function centerOfMass = computeCenterOfMass(Xv, Yv, Zv,Uvoxel,backgroundConductivity)
    
    centerOfMass = [0,0,0];

    UVoxelNew = Uvoxel-backgroundConductivity;

    UVoxelNew(isnan(UVoxelNew)) = 0;

    totalMass = sum(UVoxelNew(:));
    
    xTimesU = Xv.* UVoxelNew;
    yTimesU = Yv.* UVoxelNew;
    zTimesU = Zv.* UVoxelNew;
    
    centerOfMass(1) = sum(xTimesU(:))/totalMass;
    centerOfMass(2) = sum(yTimesU(:))/totalMass;
    centerOfMass(3) = sum(zTimesU(:))/totalMass;
end

function averagePosition = computeAveragePosition(Xv,Yv,Zv,Uvoxel,backgroundConductivity)
    
    averagePosition = [0,0,0];

    UVoxelNew = Uvoxel-backgroundConductivity;

    UVoxelNew(isnan(UVoxelNew)) = 0;

    x = Xv(UVoxelNew>1e-3);
    y = Yv(UVoxelNew>1e-3);
    z = Zv(UVoxelNew>1e-3);
    
    averagePosition(1) = mean(x(:));
    averagePosition(2) = mean(y(:));
    averagePosition(3) = mean(z(:));
end

function greitResolution = computeGreitResolution(Xv,Yv,Zv,Uvoxel,backgroundConductivity)
    
    UvoxelCorrected = computeCorrectedUvoxel(Uvoxel,backgroundConductivity);

    Vq = sum(UvoxelCorrected(not(isnan(UvoxelCorrected))));
    numberOfPixels = sum(not(isnan(Uvoxel(:))));

    greitResolution = (Vq/numberOfPixels)^(1/3);
end

function  greitShapeDeformation = computeGreitShapeDeformation(Xv,Yv,Zv,Uvoxel,backgroundConductivity)
    
    centerOfMass = computeCenterOfMass(Xv, Yv, Zv,Uvoxel,backgroundConductivity);

    UvoxelCorrected = computeCorrectedUvoxel(Uvoxel,backgroundConductivity);

    Vq = sum(UvoxelCorrected(not(isnan(UvoxelCorrected))));
    
    % C is the circle centered at centerOfMass with a volume equivalent to
    % Vq

    dx = Xv(1,2,1) - Xv(1,1,1);
    dy = Yv(2,1,1) - Yv(1,1,1);
    dz = Zv(1,1,2) - Zv(1,1,1);

    dV = dx*dy*dz;
    Vc = Vq*dV;

    R = (Vc/((4/3)*pi))^(1/3); 

    circleMask = (Xv-centerOfMass(1)).^2+(Yv-centerOfMass(2)).^2+(Zv-centerOfMass(3)).^2 <= R^2;

    pixelsInCircle = sum(UvoxelCorrected(circleMask));

    greitShapeDeformation = (Vq-pixelsInCircle)/Vq;
end

function  greitRinging = computeGreitRinging(Xv,Yv,Zv,Uvoxel,backgroundConductivity)
    
    centerOfMass = computeCenterOfMass(Xv, Yv, Zv,Uvoxel,backgroundConductivity);

    UvoxelCorrected = computeCorrectedUvoxel(Uvoxel,backgroundConductivity);

    Vq = sum(UvoxelCorrected(not(isnan(UvoxelCorrected))));
    
    % C is the circle centered at centerOfMass with a volume equivalent to
    % Vq

    dx = Xv(1,2,1) - Xv(1,1,1);
    dy = Yv(2,1,1) - Yv(1,1,1);
    dz = Zv(1,1,2) - Zv(1,1,1);

    dV = dx*dy*dz;
    Vc = Vq*dV;

    R = (Vc/((4/3)*pi))^(1/3); 

    circleMask = (Xv-centerOfMass(1)).^2+(Yv-centerOfMass(2)).^2+(Zv-centerOfMass(3)).^2 <= R^2;
    
    denominator = sum(Uvoxel(circleMask));

    notCircleMask = not(circleMask) & not(isnan(Uvoxel));

    UVoxelTemp = Uvoxel - backgroundConductivity;

    negativeMask = UVoxelTemp<0;
    
    numerator = sum(UVoxelTemp(notCircleMask & negativeMask));

    greitRinging = numerator/denominator;
end

%% Functions
function UvoxelCorrected = computeCorrectedUvoxel(Uvoxel,backgroundConductivity)
    
    % Maximum difference between anomaly conductivity and background
    maxU = max(abs(Uvoxel(:)-backgroundConductivity));

    UvoxelCorrected = Uvoxel;

    UvoxelCorrected(abs(Uvoxel-backgroundConductivity)<0.25*maxU) = 0;
    UvoxelCorrected(abs(Uvoxel-backgroundConductivity)>=0.25*maxU) = 1;
end
