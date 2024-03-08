

% Tree Construction and Node Labeling: The class constructs lineage trees by iteratively
% creating nodes representing individual cells. Each node encodes comprehensive information:
% Temporal data: Image number (time point), object number (unique cell identifier within a frame).
% Spatial data: Cell's geometric center (X, Y coordinates), and calculated depth based on the device layout.
% Biological data: Cell size (mean radius, major axis length), fluorescence intensity (e.g., GFP, RFP),
% indicating cellular states or properties.
% Lineage-specific data: Generation number, cell division cycle, lineage label, and flags for specific
% biological events (e.g., GFP rise, infection events).

% Information saved in each node is organised as below:
% 1: Image number
% 2: Object Number
% 3: GFP intensity
% 4:Radius
% 5: Length
% 6: X Location
% 7: Y Location
% (8) instantaneous (from last frame) elongation rate (zero for root)
% (9) Generation number =1 for ancestor
% (10) life history up to the point means = 1 for ancestor
% (11) GFP_rise_flag = 1 if GFP rose above threshold from
% parent to this node
% (12) infection_event_flag = 1 if infection event inferred at
% this cell-object
% (13) Depth
% (14) Chamber number
% (15) Cell division cycle
% (16) RFP intensity
% (17) Branching node
% (18) Lineage ID

global Number_Potential_Lineages
global Max_LocationXY
global  Time_Step

% Create a map where keys are column names and values are column indices
colMap = containers.Map(headers, 1:length(headers));
% identify a tree by adding a label to all the node so for example tell me
% which experiment number, tree number, chamber number does each tree have.
Max_LocationXY = [max(K(:, colMap('Location_Center_X'))), max(K(:, colMap('Location_Center_Y')))] / Pixel_to_Micrometer_Conversion;
Min_LocationXY = [min(K(:, colMap('Location_Center_X'))), min(K(:, colMap('Location_Center_Y')))] / Pixel_to_Micrometer_Conversion;

for s = 1:Number_Potential_Lineages
    r = 1;
    lineageflag = 0;

    Depth_root = Max_LocationXY(1,1) - (K(r, colMap('Location_Center_Y')/Pixel_to_Micrometer_Conversion);
    while r <= length(K)
        % recognizing cell (the ancestor => very first one) in the first frame
        if [K(r, colMap('TrackObjects_ParentImageNumber_50')),...
                K(r, colMap('TrackObjects_ParentObjectNumber_50')),...
                K(r, colMap('TrackObjects_Label_50'))] == [0,0,s]
            % Constructing the tree with the corrected column indices
            t = tree([K(r, colMap('ImageNumber')), K(r, colMap('ObjectNumber')), ...
                K(r, colMap('Intensity_MeanIntensity_GFP')),...
                K(r, colMap('AreaShape_MeanRadius')) / Pixel_to_Micrometer_Conversion, ...
                K(r, colMap('AreaShape_MajorAxisLength')) / Pixel_to_Micrometer_Conversion, ...
                K(r, colMap('Location_Center_X')) / Pixel_to_Micrometer_Conversion,...
                K(r, colMap('Location_Center_Y')) / Pixel_to_Micrometer_Conversion, ...
                NaN, 1, 1, 0, 0, Depth_root, Num_All_chambers, NaN,...
                K(r, colMap('Intensity_MeanIntensity_RFP')), 0, K(r, colMap('TrackObjects_Label_50'))]);

            lineageflag = 1;
        end
        r = r + 1;
    end

    if lineageflag == 0
        %Assigning temperary W because most trees
        % would be empty and at the end we want to filter them out
        W_temp(s) = tree;
    else

        for p = 1:size(K, 1)
            if s == K(p, colMap('TrackObjects_Label_50'))
                for u = 1:nnodes(t) % Looping through all the nodes in tree t
                    v = t.get(u); % Retrieving the contents of node u
                    w = [v(1),v(2)]; % Assigning a label to the ImageNumber and ObjectNumber of node u

                    % Check whether ImageNumber and ObjectNumber of node 'u' equal
                    % ParentImageNumber and ParentObjectNumber of cell
                    % object 'p'
                    if w == [K(p, colMap('TrackObjects_ParentImageNumber_50')), K(p, colMap('TrackObjects_ParentObjectNumber_50'))]
                        % If so: assign u as tree-parent of cell-object p
                        t = t.addnode(u, [K(p, colMap('ImageNumber')),...
                            K(p, colMap('ObjectNumber')), ...
                            K(p, colMap('Intensity_MeanIntensity_GFP')),...
                            K(p, colMap('AreaShape_MeanRadius')) / Pixel_to_Micrometer_Conversion, ...
                            K(p, colMap('AreaShape_MajorAxisLength')) / Pixel_to_Micrometer_Conversion, ...
                            K(p, colMap('Location_Center_X')) / Pixel_to_Micrometer_Conversion, ...
                            K(p, colMap('Location_Center_Y')) / Pixel_to_Micrometer_Conversion,...
                            NaN, NaN, NaN, 0, 0, 0, Num_All_chambers, NaN,...
                            K(p, colMap('Intensity_MeanIntensity_RFP')), 0, ...
                            K(p, colMap('TrackObjects_Label_50'))]);
                    end % End if for matching parent-child
                end % End for looping through nodes in tree
            end % End if for matching lineage ID
        end % End for looping through all cells
        W_temp(s) = t; % Construction of tree S complete
    end % End if for checking lineage flag
end % End for looping through potential lineages


for s = 1:Number_Potential_Lineages %loop all over lineage trees
    %adressing tree W_temp(s) in this pass through the loop
    t = W_temp(s);

    for node_index=2:nnodes(t) %Nodes number

        node_K = t.get(node_index);
        node_length = node_K(5); % cell object length of the current node
        %only for exp 4.8 change below to 1,1 and to 6, for the rest
        %should be 1,2 and 7
        % node_K(13) = Max_LocationXY(1,2) - node_K(7); %Cells Depth
        node_K(13) = Max_LocationXY(1,1) - node_K(6); %Cells Depth


        %retrieve parent of node information
        parent_of_node = t.get(t.getparent(node_index));
        parent_length = parent_of_node(5);
        parent_elongation_rate = parent_of_node(8);%elongation rate of parent node
        parent_gen_number = parent_of_node(9);% Generation number of parent node
        time_steps_since_birth_of_parent = parent_of_node(10);%number of time-steps since birth for parent node

        %next, determine generation number and timesteps_since_birth

        %check whether node corresponds to a newborn cell-object


        if length(t.getchildren(t.getparent(node_index))) == 1 %if you don't have a sister (number of daughters)
            elongation_rate = abs(node_length - parent_length)/Time_Step;  % Calculate the elongation rate in units um/min
            node_K(8) = elongation_rate;% assign elongation rate
            node_K(9) = parent_gen_number;% don't change the generation number as it hasn't been divided,
            node_K(10) = time_steps_since_birth_of_parent+1;% increase the number of time-steps since birth (i.e. lenght of cell division cycle)



        else%otherwise the node is a newborn cell-object
            node_K(8) = NaN;%no meaningful elongation rate: just born
            node_K(9) = parent_gen_number+1;% increase the generation
            node_K(10) = 1;% set the number of time-steps since birth to 1

        end
        % assign values to node
        t=t.set(node_index,node_K);
    end

    W_temp(s) = t; % All trees

end



%% Filtering out empty lineages and lineages with inappropriate length
% Filter criteria for lineage trees:
% 1) Exclude empty trees.
% 2) Exclude trees with total length less than 2µm, 
% likely due to incorrect identification of debris.
% 3) Exclude trees with total length greater than 15µm, 
% which typically represent non-viable lineages or artifacts f
% rom erroneous segmentation, such as microfluidic device walls.
% 4) Exclude trees with fewer than one node, indicative of segmentation errors.
%  This often occurs when remnants of dead cells are

% Initialize s for the current chamber
s = 1;
% Loop through each lineage in the current chamber
for k = 1:numel(W_temp)
    % Get the first node of the lineage
    V = W_temp(k).get(1);

    if ~isempty(V) && V(5) >= 2 && V(5) <= 15 && nnodes(W_temp(k)) > 2
        % If the lineage is not empty and meets the length criteria, add it to W
        W(s) = W_temp(k);
        s = s + 1; % Move to the next index within the current W cell for the same chamber
    end
end



