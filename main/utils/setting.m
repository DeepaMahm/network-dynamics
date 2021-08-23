classdef setting
    %SETTING Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        n;
    end
    
    methods
%         function data = load_setting
%             File = 'Configuration.ini';
%             I = INI('File',File);
%             I.read();
%             data = I.get('UserData'); % struct
%         end
        
        function obj = setting
            obj.n = 2;
        end
    end
end

