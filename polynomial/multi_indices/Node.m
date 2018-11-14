classdef Node < handle
    properties (SetAccess = private)
        obj
        next
        prev
    end
    
    methods
        function self = Node(obj, prev, next)
            self.obj = obj;
            self.next = next;
            self.prev = prev;
        end
        
        function self = set_next(self, next)
            if ~isa(next, 'Node')
                error("Wrong type");
            end
            self.next = next;
        end
        
        function self = set_prev(self, prev)
            if ~isa(prev, 'Node')
                error("Wrong type");
            end
            self.prev = prev;
        end
        
        function next = get_next(self)
            next = self.next;
        end
        
        function prev = get_prev(self)
            prev = self.prev;
        end
        
        function self = set_obj(self, obj)
            self.obj = obj;
        end
        
        function obj = get_obj(self)
            obj = self.obj;
        end
    end
end