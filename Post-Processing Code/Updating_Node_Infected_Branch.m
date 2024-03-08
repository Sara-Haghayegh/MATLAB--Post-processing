% In this part, we are updating the flag in element 17 of each node based on the value of element 11.


for k = 1:numel(W)
    % Get the label of the current tree
    tree_label = tree_labels(k);

    % Check if the tree is labeled as infected (label 2)
    if tree_label == 2
        % Initialize a flag to keep track if any qualifying node is found in the lineage
        lineage_has_qualifying_node = false;

        % Iterate over the nodes of the current tree to check element 11
        for i = 1:nnodes(W(k))
            node_data = W(k).get(i);

            if node_data(11) == 1
                % Set the flag in the 17th element of the current node to 1
                node_data(17) = 1;

                % Update the current node with the modified data
                W(k) = W(k).set(i, node_data);

                % Set the lineage flag to true
                lineage_has_qualifying_node = true;
            else
                % Check if the current node is a descendant of a node with element 11 equal to 1
                if lineage_has_qualifying_node
                    % Set the flag in the 17th element of the current node to 1
                    node_data(17) = 1;

                    % Update the current node with the modified data
                    W(k) = W(k).set(i, node_data);
                else
                    % Set the flag in the 17th element of the current node to 0
                    node_data(17) = 0;

                    % Update the current node with the modified data
                    W(k) = W(k).set(i, node_data);
                end
            end
        end
    end
end
