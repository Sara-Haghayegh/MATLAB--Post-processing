%% In this class, cell's state at the time of infection is investigated
% Infomation like cells size, division state and their location and time
% frame is being recorded 


% To track infection events and processed times
infectionEvents = []; % Stores frame numbers of infection events
processedInfectionTimes = []; % Tracks processed infection times to avoid redundancy

% Arrays for storing infection event-related data
Frame_Number_at_infection_event_set = [];
Length_at_infection_event_set = [];
Time_at_infection_event_set = [];
Radius_at_infection_event_set = [];
Elongation_at_infection_event_set = [];
CDCL_at_infection_event_set = []; % Cell Division Cycle Length
Generation_Number_at_infection_event_set = [];
Depth_at_infection_event_set = [];
AtDivisionState_inf = []; % At Division State for infected
Life_stage_inf = []; % Life stage for infected

% Arrays for storing non-infection event-related data at the time of infection
Length_at_Noninfection_event_set = [];
Radius_at_Noninfection_event_set = [];
Elongation_at_Noninfection_event_set = [];
CDCL_at_Noninfection_event_set = []; % Cell Division Cycle Length for non-infected
Generation_Number_at_Noninfection_event_set = [];
AtDivisionState_noninf = []; % At Division State for non-infected
Life_stage_noninf = []; % Life stage for non-infected


% Loop over all constructed lineages tree
for j = 1:numel(W)
    %% Check the infected lineages
    if tree_labels(j) == 2  % if it is infected 

        for i = 1:nnodes(W(j)) % loop over all nodes in this target_infected tree

            % retrieve node data
            infected_tree_node_data = W(j).get(i);

            if infected_tree_node_data(12) == 1  % if this cell-object corresponds to an infection event
                Frame_Number_at_infection_event_set = infected_tree_node_data(1);

                % Record the frame number of the infection event if not already done
                if ~ismember(Frame_Number_at_infection_event_set, infectionEvents)
                    infectionEvents(end + 1) = Frame_Number_at_infection_event_set;
                end

                Length_at_infection_event_set(end + 1) = infected_tree_node_data(5);
                Radius_at_infection_event_set(end + 1) = infected_tree_node_data(4);
                Generation_Number_at_infection_event_set(end + 1) = infected_tree_node_data(9);
                Frame_Number_at_infection_event_set(end + 1) = infected_tree_node_data(1);
                CDCL_at_infection_event_set(end + 1) = infected_tree_node_data(15);
                Elongation_at_infection_event_set(end + 1) = infected_tree_node_data(8);
                Depth_at_infection_event_set(end + 1) = infected_tree_node_data(13);


            end

        end
    end

    % Second: Record data for non-infected cells at infection times without redundancy
    for infectionTime = infectionEvents
        % Check if we've already processed this infection time for non-infected cells
        if ismember(infectionTime, processedInfectionTimes)
            continue; % Skip this infection time if already processed
        else
            if tree_labels(j) == 2 && infected_tree_node_data(12) == 0 && infected_tree_node_data(11) == 0 ||  tree_labels(j) == 1


                processedInfectionTimes(end + 1) = infectionTime; % Mark this infection time as processed

                for j = 1:numel(W)
                    for i = 1:nnodes(W(j))
                        node_data = W(j).get(i);

                        % Proceed only if this node's frame matches the current infection time and it's not an infected cell
                        if node_data(1) == infectionTime && node_data(12) == 0
                            Length_at_Noninfection_event_set(end + 1) = node_data(5);
                            Radius_at_Noninfection_event_set(end + 1) = node_data(4);
                            Generation_Number_at_Noninfection_event_set(end + 1) = node_data(9);
                            CDCL_at_Noninfection_event_set(end + 1) = node_data(15);
                            Elongation_at_Noninfection_event_set(end + 1) = node_data(8);
                        end
                    end
                end
            end

        end
    end
end

% Store the data for this chamber in the cell array
chamberDataAtInfection{Num_All_chambers} = struct(...
    'FrameNumberInfection', Frame_Number_at_infection_event_set, ...
    'LengthInfection', Length_at_infection_event_set, ...
    'AtDivisionStateInf', AtDivisionState_inf, ...
    'LifeStageInf', Life_stage_inf, ...
    'TimeInfection', Time_at_infection_event_set, ...
    'RadiusInfection', Radius_at_infection_event_set, ...
    'DepthInfection',Depth_at_infection_event_set,...
    'ElongationInfection', Elongation_at_infection_event_set, ...
    'CDCLInfection', CDCL_at_infection_event_set, ...
    'GenerationNumberInfection', Generation_Number_at_infection_event_set, ...
    'LengthNonInfection', Length_at_Noninfection_event_set, ...
    'AtDivisionStateNonInf', AtDivisionState_noninf, ...
    'LifeStageNonInf', Life_stage_noninf, ...
    'RadiusNonInfection', Radius_at_Noninfection_event_set, ...
    'ElongationNonInfection', Elongation_at_Noninfection_event_set, ...
    'CDCLNonInfection', CDCL_at_Noninfection_event_set, ...
    'GenerationNumberNonInfection', Generation_Number_at_Noninfection_event_set);











