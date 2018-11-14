classdef TotalOrderMultiIndexLattice < MultiIndexLattice
    properties
        levels_list  % List of lists containing references to the nodes on the same level
        init_flag = false
    end
    
    methods
        function self = TotalOrderMultiIndexLattice(dim, order)
            self = self@MultiIndexLattice(dim);
            self.levels_list = List();
            for o = 1:order+1
                self.levels_list.append( List() );
            end
        end
        
        function root = init(self)
            root = init@MultiIndexLattice(self);
            for o = 1:self.get_order()
                % To construct new level go through all the nodes in the
                % old level and add children in all the dimensions
                old_level = self.levels_list.get(o);
                cont = true;
                while cont
                    [node, cont] = old_level.next();
                    for d = 1:self.dim
                        node.add_child(d);
                    end
                end
            end
            self.init_flag = true;
        end
        
        function self = add_level(self)
            if ~self.init_flag
                error("Lattice has not been initialized yet")
            end
            order = self.get_order();
            old_level = self.levels_list.get(order+1);
            self.levels_list.append( List() );
            cont = true;
            while cont
                [node, cont] = old_level.next();
                for d = 1:self.dim
                    node.add_child(d);
                end
            end
        end
        
        function order = get_order(self)
            order = self.levels_list.size() - 1;
        end
        
        function self = append(self, node)
            append@MultiIndexLattice(self,node);
            self.levels_list.get(node.level+1).append(node);
        end
    end
end