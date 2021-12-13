classdef settings_utils
    %SETTINGS_UTILS Summary of this class goes here
    %   Detailed explanation goes here
    properties
        
        test_graph;
        run;
        
        %% bv or bv+cell
        model;
        
        %% geometry
        discretize;
        
        %% static flow calculation
        userdefined_qin; % flow input
        pbc; % pressure and flow bc
        
        %% flow input bc's
        qin;
        %% dynamic parameters
        
        % # species
        NVs; %NVs =  No. of transport species in vessel
        Vspecies;
        NCs;
        Cspecies;
        IC_Vspecies;
        IC_Cspecies;
           
        %parameters
        diffusivity;
        transporter;
        
        Vinf;
        rxn_set;
        rxn_vol;
        
        pressure_scan;
        glucose_scan;
        
        % speed up
        jpattern_set; %odeset
        jacobian_set; %odeset     
        saveop_4julia;
        
    end
    
    methods
        
        %% ---------------------------------------------------------------------------------------
        %% graph_utils settings
        %% ---------------------------------------------------------------------------------------
        
        function obj = settings_utils
            I = INI('File','input.ini');
            I.read();
            data = I.get('UserData'); % struct
            obj.test_graph = data.test_graph;
            obj.run = data.run;
            obj.model = data.model;

            obj.discretize = data.discretize;
            obj.userdefined_qin = data.userdefined_qin;
            obj.pbc = data.pbc;
            
            obj.qin = data.qin;
            
            obj.NVs = data.NVs;
            obj.Vspecies = data.Vspecies;
            obj.NCs = data.NCs;
            obj.Cspecies = data.Cspecies;
            obj.IC_Vspecies = data.IC_Vspecies;
            obj.IC_Cspecies = data.IC_Cspecies;
            
            
            obj.transporter = data.transporter;
            obj.diffusivity = data.diffusivity;
            
            obj.Vinf = data.Vinf;
            obj.rxn_set = data.rxn_set;
            obj.rxn_vol = data.rxn_vol;
            
            obj.pressure_scan = data.pressure_scan;
            obj.glucose_scan = data.glucose_scan;
            
            obj.jpattern_set = data.jpattern_set;
            obj.jacobian_set = data.jacobian_set;
            obj.saveop_4julia = data.saveop_4julia;
        end
        
    end 
end

