function h = plotSensors(sensorLocations,magnetSize)

hold on

xSensors = zeros(1,size(sensorLocations,1));
ySensors = zeros(1,size(sensorLocations,1));
zSensors = zeros(1,size(sensorLocations,1));

for sensorId = 1:size(sensorLocations,1)
    xc = sensorLocations(sensorId,1);
    yc = sensorLocations(sensorId,2);
    zc = sensorLocations(sensorId,3);

    xSensors(sensorId) = xc;
    ySensors(sensorId) = yc;
    zSensors(sensorId) = zc;

    normalVector = [xc,yc,zc]/norm([xc,yc,zc]);
    randVector = 2*(rand(1,3)-0.5);
    tangentVector1 = cross(normalVector,randVector)/norm(cross(normalVector,randVector));
    tangentVector2 = cross(normalVector,tangentVector1);

    v1 = [xc yc zc]-tangentVector1*magnetSize/2+tangentVector2*magnetSize/2;
    v2 = [xc yc zc]+tangentVector1*magnetSize/2+tangentVector2*magnetSize/2;
    v3 = [xc yc zc]+tangentVector1*magnetSize/2-tangentVector2*magnetSize/2;
    v4 = [xc yc zc]-tangentVector1*magnetSize/2-tangentVector2*magnetSize/2;

    x = [v1(1),v2(1),v3(1),v4(1)];
    y = [v1(2),v2(2),v3(2),v4(2)];
    z = [v1(3),v2(3),v3(3),v4(3)];

    % quiver3(xc,yc,zc, normalVector(1),normalVector(2), normalVector(3),'r')
    % quiver3(xc,yc,zc, tangentVector1(1), tangentVector1(2), tangentVector1(3),'g')
    % quiver3(xc,yc,zc, tangentVector2(1), tangentVector2(2), tangentVector2(3),'b')

    plot3(xc,yc,zc,'b.','MarkerSize',10);
    fill3(x,y,z,'b');
    str = strcat('S(',num2str(sensorId),')');
    text(xc,yc,zc,str,'HorizontalAlignment','left','FontSize',8);
end

% xlim([min(xSensors) max(xSensors)]);
% ylim([min(ySensors) max(ySensors)])
% zlim([min(zSensors) max(zSensors)])
end