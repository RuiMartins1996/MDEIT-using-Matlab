function  img = add_material_properties(img,conductivities)

assert(isfield(img.fwd_model,'mat_idx'),'fwd_model must contain a mat_idx field');

assert(isnumeric(conductivities));
assert(size(conductivities,1)==1 || size(conductivities,2)==1,'"conductivities" must be a vector');
assert(numel(img.fwd_model.mat_idx)==numel(conductivities),'"conductivities" must have the same number of elements has there are material regions in fwd_model');

for i = 1:numel(img.fwd_model.mat_idx)
    img.elem_data(img.fwd_model.mat_idx{i}) = conductivities(i);
end

return;

end