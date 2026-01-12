function img = mk_image_mdeit(mdl, elem_data, name)

default_name    = 'Created by mk_image_mdeit';

if nargin<3
    name = default_name;
end

if not(isstruct(mdl))
    error('Expected a struct as first argument');
end

switch mdl.type
    case 'fwd_model'
        % do nothing
    case 'image'
        if nargin == 1
            img = data_mapper(mdl);
            return
        end
        mdl = mdl.fwd_model;
    otherwise
        error('Model type must be fwd_model or image');
end


img = eidors_obj('image',name);
img.fwd_model = mdl;

if isfield(mdl,'show_slices');
    img.show_slices = mdl.show_slices;
end

img = fill_in_data(img,elem_data,no_params);

% Needed for computing magnetic field
img.fwd_solve.get_all_meas = 1;

% standard field order
img = eidors_obj('set', img);

function str = no_params
str = 'unspecified';

function img = fill_in_data(img,elem_data,params)

% If image is presented as row vector, then transpose
%  Shouldn't happen, but for matlab bugs
if size(elem_data,1) == 1 && ...
        size(elem_data,2) >  1
    elem_data = elem_data.';
end

switch size(elem_data,1)
    case 1
        sz = size(img.fwd_model.elems,1);
        elem_data_mat =    NaN*ones(sz,1);
        elem_data_mat(:) = elem_data;

        img.elem_data = elem_data_mat;
    case num_elems( img )
        img.elem_data = elem_data;
    case num_nodes( img )
        img.node_data = elem_data;
    case size(img.fwd_model.coarse2fine,2)
        img.elem_data = elem_data;
    otherwise
        error('Don''t understand number of elements.');
end
