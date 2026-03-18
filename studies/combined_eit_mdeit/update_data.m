function data_updated = update_data(data_new,file_name)

expected_fields = {'num_of_elements','num_of_measurements',...
    'rank','random_seed','sensor_radius','singular_values'};
    
    function data_old = parse_data(file_name)
        % Load the existing old data if it exists
        var = load(file_name);

        if numel(fieldnames(var)) ~= 1
            error('Expected number of structure fields to be 1');
        end
        
        field_names = fieldnames(var);

        % Fetch data_old from the only field in var
        data_old = var.(field_names{1});
    
        % Check if data_old is a structure with expected fields
        field_names = fieldnames(data_old);
        
        % Check if all expected fields are present
        all_present = all(ismember(expected_fields, field_names));

        if ~all_present
            error('data_old is missing one or more required fields.');
        end

    end

if isfile(file_name)
    
    % Fetch the old data from file
    data_old = parse_data(file_name);

    field_names = fieldnames(data_old);
    field_names(strcmp(field_names, 'singular_values')) = [];
    
    tuples_old = cell2mat( ...
        cellfun(@(f) data_old.(f)(:), field_names(:)', 'UniformOutput', false) ...
        );

    tuples_new = cell2mat( ...
        cellfun(@(f) data_new.(f)(:), field_names(:)', 'UniformOutput', false) ...
        );
    
    % Check if data_new has any new tuples missing from data_old
    is_new = ~ismember(tuples_new, tuples_old, 'rows');
    new_entry_ids = find(is_new);
    
    % Extract fields from new entries
    vals = cellfun(@(f) data_new.(f)(is_new), expected_fields, 'UniformOutput', false);
    [new_elems, new_meas, new_rank,new_rcs, new_sr] = vals{:};

    new_vals = {new_elems, new_meas,new_rank, new_rcs, new_sr};

    for k = 1:numel(field_names)
        f = field_names{k};
        data_updated.(f) = [data_old.(f); new_vals{k}];
    end
       
    num_of_old_entries = length(data_old.num_of_elements);
    num_of_new_entries = length(new_entry_ids);

    n_cols = num_of_old_entries + num_of_new_entries;
    n_rows = max(size(data_old.singular_values,1),size(data_new.singular_values,1));
    
    S = zeros(n_rows, n_cols);
    
    % Old data
    for k = 1:num_of_old_entries
        m = numel(data_old.singular_values(:, k));
        S(1:m, k) = data_old.singular_values(:, k);
    end

    % New data
    for k = 1:num_of_new_entries
        m = numel(data_new.singular_values(:, new_entry_ids(k)));
        S(1:m, num_of_old_entries + k) = data_new.singular_values(:, new_entry_ids(k));
    end

    data_updated.singular_values = S;
else
    % No existing file: all rows are new
    data_updated = data_new;
end

end