function plot_sensors(mdl)

assert(isfield(mdl,'type'),'Model must have "type" field');

switch mdl.type
    case 'image'
        mdl = mdl.fwd_model;
    case 'fwd_model'
        %do nothing
    otherwise
        error('Only "image" and "fwd_model" make sense here')
end

hold on
for m = 1:numel(mdl.sensors)

    sensor_position = mdl.sensors(m).position;
    
    xc = sensor_position(1);
    yc = sensor_position(2);
    zc = sensor_position(3);
    
    plot3(xc,yc,zc,'b.','MarkerSize',10);

    sensors_axis_1 = mdl.sensors(m).axes.axis1;
    sensors_axis_2 = mdl.sensors(m).axes.axis2;
    sensors_axis_3 = mdl.sensors(m).axes.axis3;

    % quiver3(...
    %     xc,yc,zc,...
    %     sensors_axis_1(1),sensors_axis_1(2),sensors_axis_1(3));
    % 
    % quiver3(...
    %     xc,yc,zc,...
    %     sensors_axis_2(1),sensors_axis_2(2),sensors_axis_2(3));
    % 
    % quiver3(...
    %     xc,yc,zc,...
    %     sensors_axis_3(1),sensors_axis_3(2),sensors_axis_3(3));

    str = strcat('S(',num2str(m),')');
    text(xc,yc,zc,str,'HorizontalAlignment','left','FontSize',8);

end

hold off

end