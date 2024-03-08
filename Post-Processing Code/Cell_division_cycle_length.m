

    % Loop over all constructed lineages tree
    for j = 1:numel(W)
        t = W(j); % current tree
        for node_index = 1:nnodes(t) % Nodes number of the current tree
            node_data = t.get(node_index);
            life_to_that_point = node_data(10);

            if node_data(10) == 1
                %check the time_steps_since_birth_of_parent,
                % if it is 1 means that cell divided, so check the
                % parent time_steps_since_birth_of_parent
                parent_node_index = t.getparent(node_index);

                if parent_node_index ~=0 %to prevent root (its parent index is 0)
                    parent_data = t.get (parent_node_index); % Get the parent node's data
                    % Extract the cell's life from the parent node's data
                    parent_life = parent_data(10);
                    parent_data(15) = parent_life * Time_Step; % Update the 15th element
                    % which is cell division cycle length
                    t = t.set(parent_node_index, parent_data); % Update the tree node with modified data

                end

            end
        end

        W(j) = t; % Update the lineage tree
    end