classdef MultiIndexLattice < handle
    properties
        % This list contain the MultiIndexNode objects 
        % Node 0 is the root of the lattice
        nodes_list
        dim    % dimension of the multi-indices
    end
    
    methods
        function self = MultiIndexLattice(dim)
            self.dim = dim;
            self.nodes_list = List();
        end
        
        function root = init(self)
            if self.nodes_list.size > 0
                error("Lattice already initialized")
            end
            root = MultiIndexNode(zeros(1,self.dim), self, 1);
            self.append(root);
        end
        
        function self = append(self, node)
            if ~isa(node, 'MultiIndexNode')
                error("Must be a MultiIndexNode")
            end
            self.nodes_list.append(node);
        end
        
        function size = get_size(self)
            size = self.nodes_list.size;
        end
        
        function midxs = get_midx_matrix(self)
            midxs = zeros(self.nodes_list.size, self.dim);
            cont = true;
            cntr = 0;
            while cont
                cntr = cntr + 1;
                [midx_node, cont] = self.nodes_list.next();
                midxs(cntr,:) = midx_node.idx;
            end
        end
    end
end