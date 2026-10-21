function data = GetConfigData(args)
% Reads a config. file, searches for the config. by name, which is passed
% in as an argument, and returns a struct containing the data, or an
% identifier that indicates that property is calculated as the output of a
% method specific to that configs. subclass dependent on state conditions.
%                                                    Returns data as struct
% -------------------------------------------------------------------------
% Arguments
%   1) configFile = config. file name, pass as string
%   2) configName = name of specific config. in config. file, pass as
%       string or char
% -------------------------------------------------------------------------
% Dependencies
%   #) <Dependency>
% -------------------------------------------------------------------------
% Assumptions
%   #) <Assumption>
% -------------------------------------------------------------------------
% Document Version 2.1, former versions:
%   - 2.0
%   - 1.0
% -------------------------------------------------------------------------
% MATLAB Version R2025a, also compatible with:
%   - R2024b
% -------------------------------------------------------------------------
% Developed by Alex Vance (AlexVance00 on Github)
% -------------------------------------------------------------------------

    % Allows arguments to be optional and assigned in the function call
    %   as in: FunctionTemplate(<arg_name> = <arg_val>, ...)

    % List all argument names
    arguments
        args.configFile = []; % Should NOT be vector input
        args.configName = []; % Can be vector input
    end
    arg_name_list = fieldnames(args);

    % List those argument names which are optional in 1D string array
    optional_arg_names = [];

    % Makes variables out of args' fieldnames
    for i_fieldname = 1:length(arg_name_list)
        arg_name = arg_name_list{i_fieldname};
        arg_val = args.(arg_name);

        % Input Checking
        % Checks if this argument was assigned
        if ~isempty(arg_val) 

            % Checks if assigned argument was string or char - all
            %   arguments passed should be strings or chars
            if class(arg_val) == "string" || class(arg_val) == "char"

                % Initializes assigned arguments
                eval(append(arg_name, " = arg_val;"));
            else

                % If assigned argument was not string type
                error("Passed '%s' argument must be string or char types, but was passed as %s", arg_name, class(arg_val));
            end

        % If argument was unassigned, checks if it was optional
        elseif ~ismember(arg_name, optional_arg_names)
            
            % If unassigned argument was non-optional, throws error
            error("No input for non-optional '%s' argument", arg_name);
        end
    end

    % Search MATLAB path for configFile, store matches
    configFileSearchMatches = dir(fullfile(pwd, "**", configFile));

    % Check if no matches were found
    if isempty(configFileSearchMatches)
        % If none were found, throw error- it needs to be found
        error("Config file ""%s"" not found in working directory ""%s""\n", configFile, pwd);
    end

    % If there's a match, continue. Use first match found, ignore
    % dupes
    configFilePath = fullfile(configFileSearchMatches(1).folder, configFileSearchMatches(1).name);

    splitChar = "%"; % Character denoting end of instructions section and beginning of configs.

    % Take out instructions header
    rawText = readlines(configFilePath);
    splitLineNumber = find(rawText == splitChar);
    allConfigsText = rawText(splitLineNumber + 1:end);
    % Get number of elements in configName
    % For counting them, must be strings
    if isstring(configName)
        numConfigNames = numel(configName);
    else
        numConfigNames = numel(string(configName));
    end
    data = repmat(struct, size(configName));

    % Take out everything except config. of interest, filtered by
    %   configName
    for i = 1:numConfigNames
        thisConfigName = configName(i);
        configStartLineNumber = find(allConfigsText == "\" + thisConfigName + "\");
        % Technically possible at this point for a valid string or char type
        % input to have been passed to configName that isn't in the configFile,
        % so check for that
        if isempty(configStartLineNumber)
            error("Config ""%s"" not found in ""%s""\n", thisConfigName, configFilePath);
        end
        remainingConfigsText = allConfigsText(configStartLineNumber:end);
        configEndLineNumber = find(remainingConfigsText == "\END\", 1);
        thisConfigText = remainingConfigsText(1:configEndLineNumber);
        % Remove header and footer lines
        thisConfigText = thisConfigText(2:end - 1);
        % Loop through lines
        for j = 1:length(thisConfigText)
            thisLine = thisConfigText(j);
            % Split line by spaces
            wordsOfLine = split(thisLine, " ");
            lastWord = wordsOfLine(3);
            lastValue = eval(lastWord + ";");
            % First "word" is this line's variable name
            varName = wordsOfLine(1);
            % Check if value is a string
            if ~isstring(lastValue)
                % If not, assign it to data struct
                data(i).(varName) = lastValue;
            else
                % If this line's value was a string, check for whether it was a
                % filename or "correlation"
                fileSearchMatches = dir(fullfile(pwd, "**", lastValue));
                % If any matches were folders, not files (because dir doesn't
                % know the difference), remove them
                fileSearchMatches = fileSearchMatches(~[fileSearchMatches.isdir]);
                if ~isempty(fileSearchMatches)
                    % If the string was a filename and a .xlsx file, read
                    % it in. Otherwise just set the variable value to the
                    % name of the file
                    dataFilename = fullfile(fileSearchMatches(1).folder, fileSearchMatches(1).name);
                    if ~contains(dataFilename, ".xlsx")
                        data(i).(varName) = lastValue;
                        continue;
                    end
                    cellData = readcell(dataFilename);
                    structFieldnames = cellData(1, :);
                    % Log what file was used
                    % But first check if file is a field already used
                    % in the tabulated values file
                    if ismember("file", structFieldnames)
                        error("Field ""file"" cannot be used in tabulated values at: %s\n", dataFilename);
                    end
                    data(i).(varName).file = fileSearchMatches(1).name;
                    % Check if any elements in first row are empty- this
                    % sometimes happens when excel stores empty characters
                    % in cells that used to have values but don't anymore
                    structFieldnames = structFieldnames(~ismissing(structFieldnames));
                    numColumns = numel(structFieldnames);
                    % For this varName struct, make a field for each value in
                    % structFieldnames, and assign the values under it in that
                    % column in cellData to that field
                    for k = 1:numColumns
                        % Check if value is char
                        if ischar(cellData{2, k})
                            % If so, data won't compile into correct size
                            % column of string types, it will split
                            % everything up into individual characters,
                            % with each character in its own row,
                            % regardless of how many elements are in the
                            % other columns. So make them string types.
                            % Will need to iterate thru rows
                            numRows = size(cellData, 1) - 1;
                            array = strings(numRows, 1);
                            for m = 1:numRows
                                array(m) = string(cellData(m + 1, k));
                            end
                           
                            % Now assign that array to the appropriate
                            % range in data
                            data(i).(varName).(structFieldnames{k}) = array;                            

                            % Skip to next column
                            continue;                            
                        end
                        data(i).(varName).(structFieldnames{k}) = [cellData{2:end, k}]';
                    end
    
                % Or, if the string was "correlation"
                elseif lastValue == "correlation"
                    % Set value to "flagCorrelation" flag, as a string, that
                    % will signal for the calling subclass object to
                    % initialize/overwrite this value by running its
                    % correlation method for this variable.
                    data(i).(varName) = "flagCorrelation";
                else
                    data(i).(varName) = lastValue;
                end
            end
        end
    end

    return;
end