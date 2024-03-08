%set GFP_rise_flag components 11 and 12 for nodes in target_infected trees:

global GFP_threshold
global infection_delay


for j= 1:numel(W)  %loop over all constructed lineage trees

    if tree_labels(j)==2 %target_infected trees

        %populate GFP_rise_flag components (11)
        for i = 2: nnodes(W(j))  %loop over all nodes in this target_infected tree (other than root)
            node_data=W(j).get(i); %retrieve node data
            parent_node_data=W(j).get(W(j).getparent(i));%retrieve node's tree-parent data

            if node_data(3)>GFP_threshold&& parent_node_data(3) < GFP_threshold %if node is above threshold and parent was not
                new_data = node_data;  %prepare vector of new node data
                new_data(11) = 1;   %update flag value
                W(j)=W(j).set(i,new_data);   %assign new data
            end
        end


        %populate component (12): infection_event_flag
        for i = 2: nnodes(W(j))  %loop over all nodes in this target_infected tree (other than root)
            node_data=W(j).get(i); %retrieve node data
            if node_data(11)==1  %if a GFP rise event occurred for this node
                tree_ancestor_node = i;   %initialize index to propogate back through lineage
                for k=2:infection_delay
                    if tree_ancestor_node > 1   %don't go back if at root
                        tree_ancestor_node = W(j).getparent(tree_ancestor_node);  %find parent
                    end
                end
                new_data = W(j).get(tree_ancestor_node);  %prepare vector of new node data
                new_data(12) = 1;   %update flag value
                W(j)=W(j).set(tree_ancestor_node,new_data);   %assign new data
            end
        end

    end

end

