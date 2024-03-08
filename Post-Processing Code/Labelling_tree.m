global GFP_threshold


% Initialize tree_labels for this chamber
tree_labels = zeros(1,numel(W));
% labels:
% 1: Target Non Infected
% 2: Target_Infected
% 3: Phage Producer


% Loop over all constructed lineage trees
for j = 1:numel(W)
    if isempty(W(j).get(1)) == 0 % Check if the tree is not empty
        root_node_data = W(j).get(1); % Retrieve root node data
        if root_node_data(3) > GFP_threshold % Check if root node is GFP-active
            tree_labels(j) = 3; % Define as Phage Producer
        elseif root_node_data(3) < GFP_threshold % If not GFP-active, check for target infection
            i = 1;
            while i <= nnodes(W(j))  && (isempty(tree_labels(j)) || tree_labels(j) == 0)
                node_data = W(j).get(i);
                if node_data(3) >GFP_threshold % Check if node is GFP-active
                    tree_labels(j) = 2; % Define lineage as target infected
                end
                i = i + 1;
            end
            if i == nnodes(W(j)) + 1
                tree_labels(j) = 1; % Label as uninfected target
            end
        end
    end
end
