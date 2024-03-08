
clc;
clear

%% Main Script for Image Analysis of Cell Movements and Infections
% This script is designed for the analysis of time-lapse microscopy images
% to study cell behaviors, focusing on cell movement, division,
% and infection dynamics. It utilizes lineage trees to track and analyze
% cell infection process, and behaviors over time. The tree lineage used
% in this analysis follows the structure and methodology provided in the
% documentation available at http://tinevez.github.io/matlab-tree/.
% In this work, infected cells have GFP signal and therefore code is
% adjusted based on this information.

addpath github_repo

%% Define the folder path where the Excel files are located
folderPath = '';
outputFolder =folderPath;

% Get a list of all Excel files in the folder
fileList = dir(fullfile(folderPath, '')); % Example : 'PhaseOBJ*.csv'

% Initialize arrays to store common data for all files
global Number_Potential_Lineages
global GFP_threshold
global infection_delay
global  Time_Step
global Time_threshold


%%  Defining Initial conditions of microscopy images
Time_Step = ; % Define the time interval between each images,
% Define the time length you want to process images (min)
Time_threshold = ;
% Define GFP threshold and infection delay
GFP_threshold = ;
% Delay from infection event to GFP rise, in units of time-steps
infection_delay = ; 
% Conversion factor from pixels to micrometers, based on microscopy measurements
Pixel_to_Micrometer_Conversion =;

% Initialize a cell array to store lineage data for each chamber
chamberDataAtInfection = cell(1, numel(fileList));
chamberDataOne = cell(1, numel(fileList));
chamberDataAvg = cell(1, numel(fileList));
chamberData = cell(1, numel(fileList));

%% Start the processing by defining different classes
% Loop through each Excel file/ chamber
for Num_All_chambers = 1:numel(fileList)
    % Get the current file name
    filename = fullfile(folderPath, fileList(Num_All_chambers).name);

    % Read the data from the Excel file
    K = readmatrix(filename);

    % Read the first row of the Excel file to get headers
    headers = readcell(filename, 'Range', '1:1');
    colMap = containers.Map(headers, 1:length(headers));


    % Update number of potential lineages based on the 'TrackObjects_Label_50' column
    Number_Potential_Lineages = max(K(:, colMap('TrackObjects_Label_50')));

    Buildtree
    % build trees, including elongation rate, generation, cdc length
    % produces a vector W of trees, each node has 10 componenents:
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
    % (13) Depth, calculating Depth at each time frame : Y max  -  Y current
    % node
    % (14) Chamber number: report the number of chamber in which cell is
    % loacated
    % (15) Cell division cycle
    % (16) RFP intensity
    % (17) Flagging Ancestors and Descendants with Lineage Branching

    Cell_division_cycle_length
        % calculate cell division based on the max generation number mupltiplied to
        % the time step. Update the 15th value in the tree accordingly. If cell
        % divison accures, determine the generation number for the parent and use
        % it to calculate the cell division cycle length which is :
        % Max generation number for the cell before division * Time step

    Labelling_tree
        % label trees as one of:
        % target infected (at any timepoint),
        % target non infected (otherwise)
        % labels stored in a vector tree_labels
        % labels:
        %     1 = target Not infected
        %     2 = target infected
        %     3 = Phage producer


    Infection_Flag_update
        % update infection flag
    Updating_Node_Infected_Branch
        % This part will update the flag for the infected branch. If the node is
        % infected, the node and all of its descendants' flag will be 1. If they
        % are not infected, it remains 0

    Growth

        % For all the infected lineages: All the below elements are calculated for
        %         a) each lineage
        %         b) each chamber
        %         c) each experiment
        %
        %     Elements:
        %
        %             Average and standard deviation of elongation rate
        %             Average and standard deviation of generation
        %             Average and standard deviation of Depth
        %             Avarege and standard deviation of Cell division cycle length
        %             Average and standard deviation of length
        %             Average and standard deviation of radius


      Health_at_Infection

      Distance_of_Infection
        % % % % Calculating minimum distance from phage producer cells and having
        % % % % graphs

       clear W;
       clear W_temp;
end


