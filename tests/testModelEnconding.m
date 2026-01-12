classdef testModelEnconding < matlab.unittest.TestCase
    properties
        test_parameters;
    end
    methods(TestClassSetup)
        % Shared setup for the entire test class
    end

    methods(TestMethodSetup)
        function setup(test_case)

            % Prepare workspace
            cd("C:\Users\RuiMartins\Desktop\SVD Comparison");
            addpath(genpath("functions"));
            addpath(genpath("libraries"));

            % Setup EIDORS
            eidorsFolder = setupEidors(cd);
            clc;

            run("globalParameters.m")

            % Set or create model folder
            model_folder = './models';
            if ~exist(model_folder, 'dir')
                mkdir(model_folder);
            end
            addpath(genpath("models"));

            rng(1);

            % maxsz is determined by mesh convergence study
            model_parameters.maxsz = 5;

            % This data came from Kai's Thesis
            model_parameters.isCylindrical = true; % page 163
            model_parameters.height=70; %(mm) page 163
            model_parameters.radius=40; %(mm) page 163
            model_parameters.numOfRings = 1;
            model_parameters.numOfElectrodesPerRing = 16;
            model_parameters.is2D = true;

            % what's the size of the electrodes? He mentions EEG for the tank experiment, but what for the simulated experiment?
            model_parameters.electrodeRadius = 8; %(mm) unknown from the thesis, but good estimate for Ag/AgCl EEG electrode
            assert(model_parameters.electrodeRadius*model_parameters.numOfElectrodesPerRing < 2*pi*model_parameters.radius)

            model_parameters.numOfSensors = 25; %(mm) not mentioned!
            model_parameters.sensorRadius = 70; %(mm) not mentioned!

            anomaly_conductivity = 1e-12; % (0 according to Kai's thesis, but check notes)
            anomaly_position = [20,0,35]; % pg 172
            anomaly_radius = 25/2; % 25mm of diameter, pg 172

            % This is deprecated, the .material property overrides this behaviour,
            % except for the anomaly conductivity
            model_parameters.anomaly = struct(...
                'type','spherical',...
                'conductivity',anomaly_conductivity,...
                'radius',anomaly_radius,...
                'position',anomaly_position);

            % A material defines the geometry of that material inside the domain, not
            % its properties, like conductivity.
            model_parameters.material = struct(...
                'type','cylindrical',...
                'name','plastic_cylinder',...
                'radius',anomaly_radius,...
                'position',anomaly_position);


            % Define test parameters
            test_case.test_parameters.model_parameters = model_parameters;
            test_case.test_parameters.model_folder = model_folder;
        end
    end

    methods(Test)
        % Test if built forward model is the same as the loaded forward
        % model from memory
        function test1(test_case)
            model_folder = test_case.test_parameters.model_folder;
            model_parameters = test_case.test_parameters.model_parameters;
            
            [parsed_model_parameters,fmdls] = ...
                mk_mdeit_model(model_parameters,model_folder);
            
            file_name = model_parameters_to_file_name(parsed_model_parameters,model_folder);
            var = load(file_name);

            is_same = compare_structs(var.fmdl, fmdls{1});
            
            test_case.verifyTrue(is_same);
        end

    end

end

%% FUNCTION: compareStructs
function is_same = compare_structs(a, b)
% Compare two structs recursively, ignoring row/column orientation
if ~isstruct(a) || ~isstruct(b)
    error('Inputs must be structs');
end

% Compare field names
if ~isequal(sort(fieldnames(a)), sort(fieldnames(b)))
    is_same = false;
    return;
end

% Loop over fields
flds = fieldnames(a);
for i = 1:numel(flds)
    f = flds{i};
    va = a.(f);
    vb = b.(f);

    if isstruct(va) && isstruct(vb)
        % Recursively compare structs
        if ~compare_structs(va, vb)
            is_same = false;
            return;
        end
    elseif isnumeric(va) && isnumeric(vb)
        % Compare numerics, ignoring row/column shape
        if ~isequal(va(:), vb(:))
            is_same = false;
            return;
        end
    elseif ischar(va) && ischar(vb)
        if ~strcmp(va, vb)
            is_same = false;
            return;
        end
    elseif isstring(va) && isstring(vb)
        if ~strcmp(va, vb)
            is_same = false;
            return;
        end
    elseif islogical(va) && islogical(vb)
        if ~isequal(va(:), vb(:))
            is_same = false;
            return;
        end
    else
        % Fallback strict check
        if ~isequal(va, vb)
            is_same = false;
            return;
        end
    end
end

is_same = true;
end