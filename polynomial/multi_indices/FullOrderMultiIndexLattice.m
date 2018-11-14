classdef FullOrderMultiIndexLattice < MultiIndexLattice
    properties
        levels_list  % List of lists containing references to the nodes on the same level
        max_order
        init_flag = false
    end
    
    methods
        function self = FullOrderMultiIndexLattice(dim, order)
            self = self@MultiIndexLattice(dim);
            self.max_order = order;
            self.levels_list = List();
            for o = 1:(order*self.dim+1)
                self.levels_list.append( List() );
            end
        end
        
        function root = init(self)
            root = init@MultiIndexLattice(self);
            % Not the most efficient construction but anyway...
            for o = 1:(self.get_order()*self.dim)
                % To construct new level go through all the nodes in the
                % old level and add children in all the dimensions
                old_level = self.levels_list.get(o);
                cont = true;
                while cont
                    [node, cont] = old_level.next();
                    for d = 1:self.dim
                        if node.idx(d) < self.max_order
                            node.add_child(d);
                        end
                    end
                end
            end
            self.init_flag = true;
        end
                
        function order = get_order(self)
            order = self.max_order;
        end

        function self = append(self, node)
            append@MultiIndexLattice(self,node);
            self.levels_list.get(node.level+1).append(node);
        end
    end
end