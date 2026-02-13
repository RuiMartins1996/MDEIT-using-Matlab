function plot_sensors(mdl,plot_sensor_axis,color)

assert(isfield(mdl,'type'),'Model must have "type" field');

if nargin<2
    plot_sensor_axis = false;
    color = 'default';
end

if nargin<3
    color = 'default';
end

switch mdl.type
    case 'image'
        mdl = mdl.fwd_model;
    case 'fwd_model'
        %do nothing
    otherwise
        error('Only "image" and "fwd_model" make sense here')
end

c = parse_color(color);


hold on

%Compute the bounding box of sensor positions ( for setting up axis)
max_x = max(mdl.nodes(:,1));
min_x = min(mdl.nodes(:,1));
max_y = max(mdl.nodes(:,2));
min_y = min(mdl.nodes(:,2));
max_z = max(mdl.nodes(:,3));
min_z = min(mdl.nodes(:,3));

for m = 1:numel(mdl.sensors)
    
    sensor_position = mdl.sensors(m).position;
    
    xc = sensor_position(1);
    yc = sensor_position(2);
    zc = sensor_position(3);

    max_x = max(max_x,xc);
    min_x = min(min_x,xc);
    max_y = max(max_y,yc);
    min_y = min(min_y,yc);
    max_z = max(max_z,zc);
    min_z = min(min_z,zc);
    
    plot3(xc,yc,zc,'b.','MarkerSize',10,'Color',c);

    sensors_axis_1 = mdl.sensors(m).axes.axis1;
    sensors_axis_2 = mdl.sensors(m).axes.axis2;
    sensors_axis_3 = mdl.sensors(m).axes.axis3;

    if plot_sensor_axis
        quiver3(...
            xc,yc,zc,...
            sensors_axis_1(1),sensors_axis_1(2),sensors_axis_1(3));

        quiver3(...
            xc,yc,zc,...
            sensors_axis_2(1),sensors_axis_2(2),sensors_axis_2(3));

        quiver3(...
            xc,yc,zc,...
            sensors_axis_3(1),sensors_axis_3(2),sensors_axis_3(3));
    end

    str = strcat('S(',num2str(m),')');
    text(xc,yc,zc,str,'HorizontalAlignment','left','FontSize',8);
    
    
end

hold off
axis(1.01*[min_x max_x min_y max_y min_z max_z])

end


function c = parse_color(color)

if strcmp(color,'default')
    c = 'b';
    return;
end

assert(ischar(color),'color must be a char');
assert(isscalar(color),'length of color must be 1');

allowed_colors = {'r','g','b','y','k'};
if ismember(color,allowed_colors)
    c = color;
    return;
else
    error('Color must be one of the allowed colors');
end

end