classdef dataClass < handle
    
    % A container class for holding data
    properties
        
        dataTable
        dataCount = 0;
        autoDataNameCounter = 0;
        defaultSimMin = 0.005;
        defaultSimMax = 0.7;
        
    end
    
    methods
        
        function obj = dataClass()
            
            sz = [1 4];
            varTypes = {'string','double','double','double'};
            varNames = {'Name','Data','Data Range','Simulation Range'};
            obj.dataTable = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
            
            % Set counter to zero initially. These will
            % be automatically incremented in 'appendRow'
            obj.dataCount = 0;
            obj.autoDataNameCounter = 0;
            
            % Add an initial row to the table
            obj.addData();
            
        end
        
        function obj = addData(obj,varargin)
            
            inputs = varargin{:};
            
            if isempty(inputs)
                
                    % Nothing supplied - add empty data row
                    nameVal = obj.autoDataNameCounter();
                    newName = sprintf('New data %s', nameVal);
                    
                    newDataRange = [];
                    newSimRange = [obj.defaultSimMin, obj.defaultSimMax];
                    
                    newRow = {newName, [], newDataRange, newSimRange};
                    %newUserDataRow = {'none',[0,0],[obj.defaultSimMin,obj.defaultSimMax]};
                    appendNewRow(obj,newRow);
            
            % Check length of added data
            switch length(inputs)
                case 0
                    

                    
                case 1
                    
                    % One input supplied - assume just name provided
                    newName = inputs{1};
                    if ~ischar(newName)
                        error('Single input is expected to be a data name');
                    end
                    
                    newDataRange = [];
                    newSimRange = [obj.defaultSimMin, obj.defaultSimMax];
                    
                    newRow = {newName, [], newDataRange, newSimRange};
                    %newUserDataRow = {'none',[0,0],[obj.defaultSimMin,obj.defaultSimMax]};
                    appendNewRow(obj,newRow);%,newUserDataRow);
                    
                case 2
                    
                    % Two inputs suppled - assume both name and data
                    % supplied;
                    newName = inputs{1};
                    newData = inputs{2};
                    
                    newDataX = newData(:,1);
                    newMin = newDataX(1);
                    newMax = newDataX(end);
                    
                    newDataRange = [newMin newMax];
                    newSimRange = [obj.defaultSimMin, obj.defaultSimMax];
                    
                    newRow = {newName, newData, newDataRange, newSimRange};
                    %newUserDataRow = {'data',[newMin,newMax],[newMin,newMax]};
                    appendNewRow(obj,newRow);% ,newUserDataRow);
     
                case 4
                    
                    % Four inputs = assume data and simulation ranges also
                    % supplied
                    disp('todo');
                    
                otherwise
                    
                    % Other length of inputs is not recognised
                    error('Unrecognised input into addData');
                    
            end
            
        end
        
        function displayDataObject(obj)
            
            % Display the table object. The actual obj.dataTable has the 
            % format {string, cell, double, double}, but for display we 
            % make a table that is all strings.
            
            fprintf('    Data: ------------------------------------------------------------------------------------------------------ \n\n');
            
            tab = obj.dataTable;
            

            


            % Do the custom siaplay like this*
            %tabSize = size(tab);
            
            %varTypes = {'string','string','string','string'};
            %varNames = {'Name','Data','Data Range','Simulation Range'};
            %newTable = table('Size',tabSize,'VariableTypes',varTypes,'VariableNames',varNames);
             %**Then fill thsi with strings for display**
             % e.g. 
             %  newTable{1,2} = "Data: [100 x 3]"
             %  newTable{1,3} = "[min, max]"
            
             % Create a string based table for custom display
             
             
            
            disp(tab);
 
        end
       
        
        function obj = appendNewRow(obj,newRow)
            
            tab = obj.dataTable;
            newName = newRow{1};
            if any(strcmp(newName,tab{:,1}))
                error('Duplicate data names not allowed');
            end
            
            % Carry out checks of Data type and ranges
            data = newRow{2};
            dataRange = newRow{3};
            simRange = newRow{4};
            
            if ~isempty(data)   % Data supplied
                
                if ~isnumeric(data)
                    error('Data must be a numeric array');
                end
                
                if (~isnumeric(dataRange) || (size(dataRange) ~= [1,2]) || (size(simRange) ~= [1,2]))
                    error('Data range and sim range must be [1 x 2] numeric arrays');
                end

                dataX = data(:,1);              % First column is always Q
                realDataRange = [dataX(1),dataX(end)];
                
                if dataRange(1) < realDataRange(1)
                    warning('Data range can''t be less than data min - resetting');
                    dataRange(1) = realDataRange(1);
                end
                
                if dataRange(2) > realDataRange(2)
                    warning('Data range can''t be more than data max - resetting');
                end
                
                if simRange(1) > realDataRange(1)
                    warning('Sim min range must be inside data range - resetting');
                    simRange(1) = realDataRange(1);
                end
                
                if simRange(2) < realDataRange(2)
                    warning('Sim max range must be inside data range - resetting');
                    simRange(2) = realDataRange(2);
                end
                
            else   % No data supplied
                
                dataRange = [0,0];
                
            end
            
            row = {newName, data, dataRange, simRange};
            tab = [tab ; row];
            obj.dataTable = tab;
            obj.dataCount = obj.dataCount + 1;
            obj.autoDataNameCounter = obj.autoDataNameCounter;
            
        end
        
    end
    
end
    