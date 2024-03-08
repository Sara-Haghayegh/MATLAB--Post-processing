global Time_Step;
global Time_threshold

% Initialize cell arrays to store data for each category
categories = {'Infected', 'Same_member_Not_infected',  'Non_Infected', 'Phage_Producer'};
dataByCategory_FullPeriod = struct();

% Initialize a structure to store averages
averages = struct();

for categoryIdx = 1:numel(categories)
    category = categories{categoryIdx};

    % Initialize arrays to store data for nodes in the current category
    Time = cell(1, numel(W));
    GFP = cell(1, numel(W));
    Radius_lineage = cell(1, numel(W));
    Length_lineage = cell(1, numel(W));
    Elongation_rate = cell(1, numel(W));
    Generation_data = cell(1, numel(W));
    Cell_division = cell(1, numel(W));
    Depth = cell(1, numel(W));
    RFP = cell(1, numel(W));

    % Iterate over the trees in the current chamber
    for j = 1:numel(W)
        % Initialize arrays to store data for nodes in the current tree
        numNodes = nnodes(W(j));
        Elongation = zeros(1, numNodes);
        Length = zeros(1, numNodes);
        Radius = zeros(1, numNodes);
        Cell_div = zeros(1, numNodes);
        Generation = zeros(1, numNodes);
        Depth_Experiment = zeros(1, numNodes);
        Time_Experiment = zeros(1, numNodes);
        GFP_Experiment = zeros(1, numNodes);
        RFP_Experiment = zeros(1, numNodes);

        % Iterate over the nodes of the current tree
        for i = 1:numNodes
            nodeData = W(j).get(i);

            Time_Experiment(i) = abs((nodeData(1) - First_Image_Number) * Time_Step);

            % Define a limitation for time length for better statistical analysis through all the experiments. This is optional and can be commented out.
            if Time_Experiment(i) <= Time_threshold
                % Check element 17 for category 'Infected' : Branch of
                % infected cells
                if (categoryIdx == 1 && tree_labels(j) == 2 && nodeData(17) == 1) % Infected
                    Elongation(i) = nodeData(8);
                    Length(i) = nodeData(5);
                    Radius(i) = nodeData(4);
                    Cell_div(i) = nodeData(15);
                    Generation(i) = nodeData(9);
                    Depth_Experiment(i) = nodeData(13);
                    GFP_Experiment(i) = nodeData(3);
                    RFP_Experiment(i) = nodeData(16);
                elseif (categoryIdx == 2 && tree_labels(j) == 2 && nodeData(17) ~= 1) % Same_member_Not_infected
                    Elongation(i) = nodeData(8);
                    Length(i) = nodeData(5);
                    Radius(i) = nodeData(4);
                    Cell_div(i) = nodeData(15);
                    Generation(i) = nodeData(9);
                    Depth_Experiment(i) = nodeData(13);
                    GFP_Experiment(i) = nodeData(3);
                    RFP_Experiment(i) = nodeData(16);
                elseif (categoryIdx == 4 && tree_labels(j) == 3) % Phage_Producer
                    Elongation(i) = nodeData(8);
                    Length(i) = nodeData(5);
                    Radius(i) = nodeData(4);
                    Cell_div(i) = nodeData(15);
                    Generation(i) = nodeData(9);
                    Depth_Experiment(i) = nodeData(13);
                    GFP_Experiment(i) = nodeData(3);
                    RFP_Experiment(i) = nodeData(16);
                elseif (categoryIdx == 3 && tree_labels(j) == 1) % Non_Infected
                    Elongation(i) = nodeData(8);
                    Length(i) = nodeData(5);
                    Radius(i) = nodeData(4);
                    Cell_div(i) = nodeData(15);
                    Generation(i) = nodeData(9);
                    Depth_Experiment(i) = nodeData(13);
                    GFP_Experiment(i) = nodeData(3);
                    RFP_Experiment(i) = nodeData(16);
                end
            end
        end

        % Check if the lineage has any non-empty data
        if any(Elongation) || any(Length) || any(Radius) || any(Cell_div) || any(Generation) || any(Depth_Experiment) || any(Time_Experiment) || any(GFP_Experiment) || any(RFP_Experiment)
            % Store data for the current tree
            Elongation_rate{j} = Elongation;
            Generation_data{j} = Generation;
            Length_lineage{j} = Length;
            Radius_lineage{j} = Radius;
            Cell_division{j} = Cell_div;
            Depth{j} = Depth_Experiment;
            Time{j} = Time_Experiment;
            GFP{j} = GFP_Experiment;
            RFP{j} = RFP_Experiment;
        end
    end

    % Store lineage data for this category
    dataByCategory_FullPeriod.(category) = struct( ...
        'Elongation_rate', Elongation_rate', ...
        'Generation_data', Generation_data', ...
        'Length_lineage', Length_lineage', ...
        'Radius_lineage', Radius_lineage', ...
        'Cell_division', Cell_division', ...
        'Depth', Depth', ...
        'Time', Time', ...
        'GFP', GFP', ...
        'RFP', RFP' ...
        );

    % Calculate the average for each chamber
    averages.(category) = calculateAverages(dataByCategory_FullPeriod.(category));

end

% Store the data for each category in chamberData
chamberDataOne{Num_All_chambers} = dataByCategory_FullPeriod;
chamberDataAvg{Num_All_chambers} = averages;

%% Function to calculate averages
function avgData = calculateAverages(data)
avgData = struct();

% Iterate over the fields in the data structure
fieldNames = fieldnames(data);
for i = 1:numel(fieldNames)
    fieldName = fieldNames{i};
    fieldData = data.(fieldName);

    % Check if the data is a cell array
    if iscell(fieldData)
        % Calculate the average for each chamber
        avgData.(fieldName) = cellfun(@(x) nanmean(cell2mat(x), 2), fieldData, 'UniformOutput', false);
    else
        % Calculate the average for numeric arrays
        avgData.(fieldName) = nanmean(fieldData, 2);
    end
end
end
