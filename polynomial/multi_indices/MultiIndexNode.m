classdef MultiIndexNode < handle
    properties (SetAccess = private)
        idx            % array representing the multi-index
        level          % level to which the multi-index belongs (level = sum(idx))
        dim            % dimension of the multi-index (dim = size(idx,2))
        lattice        % overall lattice object
        lattice_idx    % index position of the node in lattice.nodes_list
        parents_list   % list parent nodes
        children_list  % list of child nodes
    end
    
    methods
        function self = MultiIndexNode(idx, lattice, lattice_idx)
            if size(idx,2) ~= lattice.dim
                error("Dimension mismatch between index and lattice")
            end
            % Set properties
            self.idx = idx;
            self.level = sum(idx);
            self.dim = size(idx,2);
            self.lattice = lattice;
            self.lattice_idx = lattice_idx;
            % Set list of children to empty
            self.children_list = List();
            for d = 1:self.dim
                self.children_list.append(false);
            end
            % Set list of parents to an empty list of size dim (if level <
            % dim, some of the parents will stay empty)
            self.parents_list = List();
            for d = 1:self.dim
                self.parents_list.append(false);
            end        
        end
        
        function self = add_child(self, d)
            % Add a child with increased multi-index in dimension d
            if d > size(self.idx,2)
                error("Invalid increasing dimension")
            end
            if self.children_list.get(d) == false  % If the child does not exist already
                new_idx = self.idx;
                new_idx(d) = new_idx(d) + 1;
                child = MultiIndexNode(new_idx, self.lattice, ...
                    self.lattice.get_size()+1);
                % Append to the lattice.nodes_list 
                self.lattice.append(child);
                % Add the child to the list of children
                self.children_list.set(d, child);
                % Add self to the list of parents of the child
                child.parents_list.set(d, self);
                % For all parents p of self, connect the d-th child (uncle)
                % to child.
                cont = true;
                pn = 0;
                while cont
                    % The number pn indicates the placement in the
                    % child.parents_list and uncle.children_list
                    [p, cont] = self.parents_list.next();
                    pn = pn + 1;  
                    if ~isequal(p, false) % Empty parent
                        uncle = p.children_list.get(d);
                        uncle.children_list.set(pn, child);
                        child.parents_list.set(pn, uncle);
                    end
                end
            end
        end
        
        function idx = get_idx(self)
            idx = self.get_idx;
        end
        
        function l = get_level(self)
            l = sum(self.idx);
        end
        
        function ps = get_parents(self)
            ps = self.parents_list;
        end
        
        function cs = get_children(self)
            cs = self.children_list;
        end
        
    end
end