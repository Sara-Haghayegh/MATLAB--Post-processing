%% In this class, distance of host cells either infected or non infected from host is calculated
% Three types of cells object was recognized: 
% 1) Non Infected cells: 
    % a) From the same lineage as infected cells 
    % b) hosts with no relation to infected cells and non infected  
% 2) Infected cells
% 3) Phage Producer cells


%determine distances from infection events to producer cells

Mindistance_to_infection_event_set = [];
Sumdistance_to_infection_event_set = [];
num_producers_to_infection_event_set = [];
infection_Frame_Number = [];

Mindistance_to_noninfection_event_set = [];
Sumdistance_to_noninfection_event_set = [];
num_producers_to_noninfection_event_set = [];

Mindistance_to_uninfected_lineage_event_set = [];
Sumdistance_to_uninfected_lineage_event_set = [];
num_producers_to_uninfected_lineage_event_set = [];



% Initialize counters
mindistance_inf = 1;
mindistance_noninf = 1;
mindistance_healthy = 1;



for j = 1:numel(W)  % loop over all constructed lineage trees

    if tree_labels(j) == 2 || tree_labels(j) == 1 % target_infected trees

        for i = 1:nnodes(W(j))  % loop over all nodes in this target_infected tree

            recipient_tree_node_data = W(j).get(i); % retrieve node data
            recipient_framenumber = recipient_tree_node_data(1);
            recipient_location_X = recipient_tree_node_data(6);
            recipient_location_Y = recipient_tree_node_data(7);
            minDistance = Inf;
            sum_exp_neg_distance = 0;
            num_producer_neighbors = 0;

            for k = 1:numel(W)  % loop over all constructed lineage trees

                if tree_labels(k) == 3 % Phage producer trees

                    for kk = 1:nnodes(W(k))  % loop over all nodes in this producer tree

                        producer_node_data = W(k).get(kk); % retrieve node data

                        if producer_node_data(1) == recipient_framenumber  % if producer cell-object shares framenumber with infection event
                            producer_location_X = producer_node_data(6);
                            producer_location_Y = producer_node_data(7);
                            distance = sqrt((recipient_location_X - producer_location_X)^2 + (recipient_location_Y - producer_location_Y)^2);
                            minDistance = min(minDistance, distance);
                            sum_exp_neg_distance = sum_exp_neg_distance + exp(-distance);
                            num_producer_neighbors = num_producer_neighbors + 1;

                        end
                    end
                end
            end
            if recipient_tree_node_data(12) == 1% if this cell-object corresponds to an infection event
                infection_Frame_Number (end+1) = recipient_framenumber;
                % Store data in the arrays

                Mindistance_to_infection_event_set(end+1) = minDistance;
                Sumdistance_to_infection_event_set(end+1) = sum_exp_neg_distance;
                num_producers_to_infection_event_set(end+1) = num_producer_neighbors;

                % Increment counters
                mindistance_inf = mindistance_inf + 1;
            end


            match = false; % Define a flag to check the frame number. update it if it has similar frame number
            for k = 1: numel(infection_Frame_Number)
                if recipient_tree_node_data(1) == infection_Frame_Number(k)
                    match = true;
                    break;
                end
            end

            if match
                if  recipient_tree_node_data(12) == 0 && recipient_tree_node_data(11) == 0  % if this cell-object doesn't correspond to an infection event

                    % Store data in the arrays
                    Mindistance_to_noninfection_event_set(end+1) = minDistance;
                    Sumdistance_to_noninfection_event_set(end+1) = sum_exp_neg_distance;
                    num_producers_to_noninfection_event_set(end+1) = num_producer_neighbors;

                    % Increment counters
                    mindistance_noninf = mindistance_noninf + 1;
                end

                if tree_labels(j)==1 && recipient_tree_node_data(1) == recipient_framenumber %if this cell-object is in a non 'infected' lineage

                    Mindistance_to_uninfected_lineage_event_set(end+1) = minDistance;
                    Sumdistance_to_uninfected_lineage_event_set(end+1) = sum_exp_neg_distance;
                    num_producers_to_uninfected_lineage_event_set(end+1) = num_producer_neighbors;
                    mindistance_healthy = mindistance_healthy + 1;
                end
            end
        end

    end
end




% Store lineage data for this chamber in chamberData
chamberData{Num_All_chambers} = struct( ...
    'Mindistance_to_infection_event_set', Mindistance_to_infection_event_set, ...
    'Sumdistance_to_infection_event_set', Sumdistance_to_infection_event_set, ...
    'num_producers_to_infection_event_set', num_producers_to_infection_event_set, ...
    'Mindistance_to_noninfection_event_set', Mindistance_to_noninfection_event_set, ...
    'Sumdistance_to_noninfection_event_set', Sumdistance_to_noninfection_event_set, ...
    'num_producers_to_noninfection_event_set', num_producers_to_noninfection_event_set, ...
    'Mindistance_to_uninfected_lineage_event_set', Mindistance_to_uninfected_lineage_event_set, ...
    'Sumdistance_to_uninfected_lineage_event_set', Sumdistance_to_uninfected_lineage_event_set, ...
    'num_producers_to_uninfected_lineage_event_set', num_producers_to_uninfected_lineage_event_set ...
    );

% Clear the arrays for the next chamber
Mindistance_to_infection_event_set = [];
Sumdistance_to_infection_event_set = [];
num_producers_to_infection_event_set = [];
Mindistance_to_noninfection_event_set = [];
Sumdistance_to_noninfection_event_set = [];
num_producers_to_noninfection_event_set = [];
Mindistance_to_uninfected_lineage_event_set = [];
Sumdistance_to_uninfected_lineage_event_set = [];
num_producers_to_uninfected_lineage_event_set = [];





