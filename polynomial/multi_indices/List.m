classdef List < handle
    properties (SetAccess = private)
        head
        bottom
        size
        iter_pntr
    end
    
    methods 
        function self = List()
            self.head = false;
            self.bottom = false;
            self.size = 0;
            self.iter_pntr = false;
        end
        
        function self = append(self, obj)
            new_node = Node(obj, self.bottom, false);
            if isequal(self.bottom, false)
                self.head = new_node;
                self.bottom = new_node;
            else
                self.bottom.set_next( new_node );
                self.bottom = new_node;
            end
            self.size = self.size + 1;
        end
        
        function obj = get(self, i)
            if i >  self.size
                error("Index out of bound")
            end
            n = self.head;
            for k = 2:i
                n = n.get_next();
            end
            obj = n.get_obj();
        end
        
        function set(self, i, obj)
            if i >  self.size
                error("Index out of bound")
            end
            n = self.head;
            for k = 2:i
                n = n.get_next();
            end
            n.set_obj(obj);
        end
        
        function [obj, cont] = next(self)
            if isequal(self.iter_pntr, false)
                self.iter_pntr = self.head;
            end
            if isequal(self.iter_pntr, false)
                obj = false;
                cont = false;
            else
                nxt = self.iter_pntr.get_next();
                cont = ~isequal(nxt, false);
                obj = self.iter_pntr.get_obj();
                self.iter_pntr = nxt;
            end
        end
        
        function reset(self)
            self.iter_pntr = false;
        end
    end
end