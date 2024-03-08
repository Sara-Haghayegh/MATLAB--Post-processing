%% In this class, cell's state at the time of infection is investigated
% Infomation like cells size, division state and their location and time
% frame is being recorded 


% Initialize empty arrays for Length, Radius, Elongation, and CDCL
Frame_Number_at_infection_event_set = [];
Length_at_infection_event_set = [];
Length_at_Noninfection_event_set = [];
Length_Healthy_Cell_at_infection_event_set = [];
Length_Phage_producer_at_infection_event_set = [];
Time_at_infection_event_set = [];
AtDivisionState_inf = [];
Life_stage_inf = [];
AtDivisionState_noninf = [];
Life_stage_noninf = [];
AtDivisionState_healthy = [];
Life_stage_healthy = [];
Depth_at_infection_event_set =[];


Radius_at_infection_event_set = [];
Radius_at_Noninfection_event_set = [];
Radius_Healthy_Cell_at_infection_event_set = [];
Radius_Phage_producer_at_infection_event_set = [];

Elongation_at_infection_event_set = [];
Elongation_at_Noninfection_event_set = [];
Elongation_Healthy_Cell_at_infection_event_set = [];
Elongation_Phage_producer_Cell_at_infection_event_set = [];

CDCL_at_infection_event_set = [];
CDCL_at_Noninfection_event_set = [];
CDCL_Healthy_Cell_at_infection_event_set = [];
CDCL_Phage_producer_Cell_at_infection_event_set = [];

Generation_Number_at_infection_event_set = [];
Generation_Number_at_Noninfection_event_set = [];
Generation_Number_Healthy_Cell_at_infection_event_set = [];
Generation_Number_Phage_producer_Cell_at_infection_event_set = [];



% Loop over all constructed lineages tree
for j = 1:numel(W)
    %% Check the infected lineages health
    if tree_labels(j) == 2 % if it is infected

        for i = 1:nnodes(W(j)) % loop over all nodes in this target_infected tree

            % retrieve node data
            infected_tree_node_data = W(j).get(i);

            if infected_tree_node_data(12) == 1  % if this cell-object corresponds to an infection event
                Frame_Number_at_infection_event_set = infected_tree_node_data(1);

                Length_at_infection_event_set(end + 1) = infected_tree_node_data(5);
                Radius_at_infection_event_set(end + 1) = infected_tree_node_data(4);
                Generation_Number_at_infection_event_set(end + 1) = infected_tree_node_data(9);
                Frame_Number_at_infection_event_set(end + 1) = infected_tree_node_data(1);
                CDCL_at_infection_event_set(end + 1) = infected_tree_node_data(15);
                Elongation_at_infection_event_set(end + 1) = infected_tree_node_data(8);
                Depth_at_infection_event_set(end + 1) = infected_tree_node_data(13);
                Time_at_infection_event_set(end + 1) = (infected_tree_node_data(1)-First_Image_Number(Num_All_chambers));


                %check the life stage of node, if it is newborn
                %node(10) should be 1 otherwise it can have other
                %values
                Life_stage_inf(end+1) = infected_tree_node_data(10);
                % Check if the node has daughters
                daughters = W(j).getchildren(i);

                if length(daughters) > 0
                    AtDivisionState_inf(end + 1) = 1;
                else
                    AtDivisionState_inf(end + 1) = 0;
                end

            end



            if infected_tree_node_data(12) == 0 && infected_tree_node_data(11) == 0  % if this cell-object doesn't correspond to an infection event
                % Initialize match as false
                match = false;

                % Check if the frame number matches any infection event frame number
                for k = 1:numel(Frame_Number_at_infection_event_set)
                    if infected_tree_node_data(1) == Frame_Number_at_infection_event_set(k)
                        % Set match to true and break the loop
                        match = true;
                        break;
                    end
                end

                % If there's a match, add the length to the non-infection event set
                if match
                    Length_at_Noninfection_event_set(end + 1) = infected_tree_node_data(5);
                    Radius_at_Noninfection_event_set(end + 1) = infected_tree_node_data(4);
                    Generation_Number_at_Noninfection_event_set(end + 1) = infected_tree_node_data(9);
                    CDCL_at_Noninfection_event_set(end + 1) = infected_tree_node_data(15);
                    Elongation_at_Noninfection_event_set(end + 1) = infected_tree_node_data(8);

                    Life_stage_noninf(end+1) = infected_tree_node_data(10);
                    % Check if the node has daughters
                    daughters = W(j).getchildren(i);

                    if length(daughters) > 0
                        AtDivisionState_noninf(end + 1) = 1;
                    else
                        AtDivisionState_noninf(end + 1) = 0;
                    end
                end
            end
        end
    end

    %% Now check the states of cells that didn't get infected throughout the experiment
    if tree_labels(j) == 1 % if it is Not infected

        for i = 1:nnodes(W(j)) % loop over all nodes in this target_NonInfected tree

            % retrieve node data
            NonInfected_tree_node_data = W(j).get(i);
            Non_Infected_framenumber = NonInfected_tree_node_data(1);

            for k = 1:numel(Frame_Number_at_infection_event_set)
                if Non_Infected_framenumber == Frame_Number_at_infection_event_set(k)   % if this cell-object is in the same frame as the infected object
                    Length_Healthy_Cell_at_infection_event_set(end + 1) = NonInfected_tree_node_data(5);
                    Radius_Healthy_Cell_at_infection_event_set(end + 1) = NonInfected_tree_node_data(4);
                    Generation_Number_Healthy_Cell_at_infection_event_set(end + 1) = NonInfected_tree_node_data(9);
                    CDCL_Healthy_Cell_at_infection_event_set(end + 1) = NonInfected_tree_node_data(15);
                    Elongation_Healthy_Cell_at_infection_event_set(end + 1) = NonInfected_tree_node_data(8);
                    Life_stage_healthy(end+1) = infected_tree_node_data(10);
                    % Check if the node has daughters
                    daughters = W(j).getchildren(i);

                    if length(daughters) > 0
                        AtDivisionState_healthy(end + 1) = 1;
                    else
                        AtDivisionState_healthy(end + 1) = 0;
                    end

                end
            end
        end
    end

    %% Now check the states of phage producer cells at the time of infection
    if tree_labels(j) == 3 % if it is Phage Producer Cell

        for i = 1:nnodes(W(j)) % loop over all nodes in this target_PhageProducer tree

            % retrieve node data
            Phage_Producer_tree_node_data = W(j).get(i);

            for k = 1:numel(Frame_Number_at_infection_event_set)
                if Phage_Producer_tree_node_data(1) == Frame_Number_at_infection_event_set(k)   % if this cell-object is in the same frame as the infected object
                    Length_Phage_producer_at_infection_event_set(end + 1) = Phage_Producer_tree_node_data(5);
                    Radius_Phage_producer_at_infection_event_set(end + 1) = Phage_Producer_tree_node_data(4);
                    Generation_Number_Phage_producer_Cell_at_infection_event_set(end + 1) = Phage_Producer_tree_node_data(9);
                    CDCL_Phage_producer_Cell_at_infection_event_set(end + 1) = Phage_Producer_tree_node_data(15);
                    Elongation_Phage_producer_Cell_at_infection_event_set(end + 1) = Phage_Producer_tree_node_data(8);
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
    'GenerationNumberNonInfection', Generation_Number_at_Noninfection_event_set, ...
    'LengthHealthy', Length_Healthy_Cell_at_infection_event_set, ...
    'AtDivisionStateHealthy', AtDivisionState_healthy, ...
    'LifeStageHealthy', Life_stage_healthy, ...
    'RadiusHealthy', Radius_Healthy_Cell_at_infection_event_set, ...
    'ElongationHealthy', Elongation_Healthy_Cell_at_infection_event_set, ...
    'CDCLHealthy', CDCL_at_Noninfection_event_set, ...
    'GenerationNumberHealthy', Generation_Number_Healthy_Cell_at_infection_event_set, ...
    'LengthPhage', Length_Phage_producer_at_infection_event_set, ...
    'RadiusPhage', Radius_Phage_producer_at_infection_event_set, ...
    'ElongationPhage', Elongation_Phage_producer_Cell_at_infection_event_set, ...
    'CDCLPhage', CDCL_Phage_producer_Cell_at_infection_event_set, ...
    'GenerationNumberPhage', Generation_Number_Phage_producer_Cell_at_infection_event_set ...
    );

