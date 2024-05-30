# Perception Tools Dev Log

| Pengfei LIU |
| ----------- |

## MDF signal conversion

Before using the Perception tool, it is necessary to convert the raw signals recorded in MF4 format by CANape into .mat format, which can be directly used by MATLAB. This conversion is achieved through the `convertLog.m` script, which utilizes the `Read_MDF_DM.m`, `readMDF_DM_MD.m`, and Toolbox functions.

For the SWEET 200 signals, there are no issues. When converting SWEET 420 signals using this script, there are cases where signal names longer than 63 characters are truncated. This issue is from MATLAB's limitation on variable name lengths, with some SWEET 420 signal names being quite long.

To address this issue, I referred to the algorithm in `Read_MDF_DM_AC.m` to shorten the signal names. The approach is as follows:

- Checking if the signal name contains a dot `.`:
  - If the signal name contains a dot, it is split into different parts `(strsplit(curNames{i}{j}, '.'))`. Any possible square brackets in each part are removed before recombining them. Typically, the last two parts are combined, especially if their combined length is less than 62 characters. If this combined length is still too long, only the last part is used as the new signal name.

- Handling brackets:
  - Replace left brackets `[` with an underscore `_`.
  - Remove right brackets `]`.

- Handling non-dotted and overly long names:
  - If the signal name does not contain a dot and exceeds 63 characters, the last 63 characters of the signal name are directly used.


The algorithm in `Read_MDF_DM_AC.m`:

```matlab
for i=1:length(curNames(:,1))
    timeVectors{i, 1} = ['t' num2str(i)];
    try timeVectors{i, 2} = read(mdfObj, i, mdfObj.ChannelNames{i}{curSizes(i)}, 'OutputFormat','Vector');
        if isempty(timeVectors{i,2})
            timeVectorsErrors(i, 1) = 1;
        end
    catch
        timeVectorsErrors(i, 1) = 1;
    end
    for j=1:(curSizes(i) - 1) 
        % ici on check chaque signal du groupe de signaux : 
        
        if length(curNames{i}{j})>63 %Si le nom est trop fat, on enlève les préfix
            if contains(curNames{i}{j},'.')
                temp = strsplit(curNames{i}{j},'.');
                temp2 = string(temp);
                temp2 = strrep(strrep(temp2,']',''),'[','_');

                if (length(char(temp2(:,end))) + length(char(temp2(:,end-1))))<62
                    curNames{i}{j} = [char(temp2(:,end-1)) '_' char(temp2(:,end))];
                else
                    curNames{i}{j} = char(temp(:,end));
                end
            else % si il n'y a pas de préfix, on guillotine, gardant la fin pour distinguer chaque signaux correctement
                curNames{i}{j}(length(curNames{i}{j})-63:end);
            end
        end
        if contains(curNames{i}{j},'[')
            curNames{i}{j} = strrep(curNames{i}{j},'[','_');

        end
        if contains(curNames{i}{j},']')
                       curNames{i}{j} = strrep(curNames{i}{j},']','');
        end
            if isempty(listVars) || any(contains(listVars,strrep(curNames{i}{j},'.','_'))) % if variable list not define OR if current signal found in specified variable list
            if (Interpolation)
                dataTable{last_index, 1} = curNames{i}{j};
            else
                dataTable{last_index, 1} = [curNames{i}{j} '_' timeVectors{i, 1}];
            end
            dataTable{last_index, 3} = timeVectors{i,2};
            try dataTable{last_index, 2} = read(mdfObj, i, mdfObj.ChannelNames{i}{j}, 'OutputFormat','Vector');
            catch
                dataTableErrors(last_index, 1) = 1;
            end
            last_index = last_index + 1;
        end
    end
end
```

But the rest of this version of the code is overly simplified, particularly lacking error handling and feedback. Thus I combined the original `readMDF_DM_MD.m` and `Read_MDF_DM_AC.m` to update the code, making improvements and extensions:

- Supports fuzzy matching to filter based on the start, end, or specific patterns of signal names. This is achieved by adding a new `wildcard` function, which can handle variable lists marked with wildcards like `*`, as well as functions like `listStatsWith`, `listEndsWith`, and `listStartAndEndsWith`.
- Integrates a `waitbar`.
- Significantly improves error handling, e.g., capturing and displaying error messages when MDF file reading fails, then closing the waitbar and exiting the function.
- Adjusts details in handling overly long signal names, particularly managing character lengths more accurately when dealing with dot-separated names.
- Customizes shortening rules for signal names further through the use of regexp and contains functions, dynamically adjusting signal names based on specified name patterns in the listNames parameter.
- Extends data interpolation and time handling:
  - In interpolation mode, the code generates a new time vector to ensure all signals are interpolated at the same time points, which is crucial for subsequent data analysis.
  - Adds handling for POSIX time, including the initial timestamp's corresponding POSIX time in the final data table, which is beneficial for further analysis of time series data.
- Enhances code maintainability and extensibility:
  - The function structure is clearer, and the introduction of error handling and user feedback mechanisms makes the code more robust and easier to maintain.
  - The new fuzzy matching feature and progress bar mechanism allow for easy modifications or enhancements as needed.

**Source code of new `Read_MDF_DM.m`:**

```matlab
function [finalDatas,  FieldMatrix] = Read_MDF_DM(curPath, Interpolation, SampleTime, SetStartTime0,listVars,listNames)
if nargin < 6 % listNames not defined
    listNames = {};
end
if nargin < 5 % listVars not defined
    listVars = {};
end

listStatsWith = cellfun(@(x) x(1:end-1),listVars(endsWith(listVars,'*')),'UniformOutput',false);
listEndsWith  = cellfun(@(x) x(2:end),listVars(startsWith(listVars,'*')),'UniformOutput',false);
listStartAndEndsWith = cellfun(@(x) x(2:end-1),listVars(startsWith(listVars,'*')&endsWith(listVars,'*')),'UniformOutput',false);
listMatch     = listVars(~contains(listVars,'*'));


FieldMatrix = {};
finalDatas = {};
% Read_MDF

% [FileName,PathName,~] = uigetfile('*.mdf;*.mf4','MultiSelect','on');
warning('off','all');
curBarH = waitbar(0,'Parsing MDF...');
curBarHb=findobj(curBarH,'Type','figure');
curBarHt = get(get(curBarHb,'currentaxes'),'title');
try
    mdfObj = mdf(curPath);
catch ME
    warning(ME.message);
    fprintf('\t --> Ignored.\n');
    close(curBarH);
    return;
end

% Parse all available rasters
curNames = get(mdfObj, 'ChannelNames');
curSizes = cellfun(@length, curNames);
timeVectors = cell(length(curNames(:,1)), 2);
timeVectorsErrors = zeros(length(curNames(:,1)), 1);
dataTable = cell(sum(curSizes), 3);
dataTableErrors = zeros(sum(curSizes), 1);
last_index = 1;
waitbar(0.1,curBarH);

%% Modified for signal name more than 63 char  
for i=1:length(curNames(:,1))
    timeVectors{i, 1} = ['t' num2str(i)];
    try timeVectors{i, 2} = read(mdfObj, i, mdfObj.ChannelNames{i}{curSizes(i)}, 'OutputFormat','Vector');
        if isempty(timeVectors{i,2})
            timeVectorsErrors(i, 1) = 1;
        end
    catch
        timeVectorsErrors(i, 1) = 1;
    end
 
    for j=1:(curSizes(i) - 1)
        curSignalName = curNames{i}{j};
        if length(curSignalName) > 63 % Apply the method of 'Read_MDF_DM_AC.m' 
            if contains(curSignalName,'.')
                parts = strsplit(curSignalName,'.');
                parts = strrep(strrep(parts,']',''),'[','_');
                if (length(char(parts(end))) + length(char(parts(end-1))) < 62)
                    curSignalName = [char(parts(end-1)) '_' char(parts(end))];
                else
                    curSignalName = char(parts(end));
                end
            else
                curSignalName = curSignalName(end-62:end);
            end
        end
        
        curNames{i}{j} = strrep(curSignalName, '[', '_');  
        curNames{i}{j} = strrep(curNames{i}{j}, ']', '');  

        if isempty(listVars) ||  wildcard(strrep(curNames{i}{j},'.','_'),listStatsWith,listEndsWith,listStartAndEndsWith,listMatch)
            if (Interpolation)
                dataTable{last_index, 1} = curNames{i}{j};
            else
                dataTable{last_index, 1} = [curNames{i}{j} '_' timeVectors{i, 1}];
            end
            dataTable{last_index, 3} = timeVectors{i,2};
            try dataTable{last_index, 2} = read(mdfObj, i, mdfObj.ChannelNames{i}{j}, 'OutputFormat','Vector');
            catch
                dataTableErrors(last_index, 1) = 1;
            end
            last_index = last_index + 1;
        end
    end
    waitbar(0.1+i/length(curNames(:,1))*0.8,curBarH);
end

dataTable = dataTable(~dataTableErrors, :);
dataTable = dataTable(~cellfun(@isempty, dataTable(:, 1)), :);
timeVectors = timeVectors(~timeVectorsErrors, :);
timeVectors = timeVectors(~cellfun(@isempty, timeVectors(:, 1)), :);

curTimes2Print = cellfun(@(x,y) [x,'_Length_of_t=',num2str(length(y))],timeVectors(:,1),timeVectors(:,2),'UniformOutput',false);

if (~Interpolation)
    [ChosenTimeIdx,~] = listdlg('ListString',curTimes2Print,'Name','Temps','ListSize',[300 400]);
    t=timeVectors{ChosenTimeIdx,2}(:,1);
    dataTable = dataTable(:,1:2);
    [~,idx] = sort(upper(dataTable(:,1)));
    finalDatas = [dataTable(idx, :) ; timeVectors ; {'t', t}];
else
    curBarHt.String = 'Interpolation...';
    ValidIdx = cellfun(@(x, y) ~isempty(x) && ~isempty(y) && length(x) == length(y), dataTable(:,2) , dataTable(:,3));
    MaxT = max(cellfun(@(x) x(end), timeVectors(:,2)));
    MinT = min(cellfun(@(x) x(1), timeVectors(:,2)));
    NewTime = (MinT:SampleTime:MaxT)';
    dataTable = dataTable(ValidIdx, :);
    %     dataTable(:,3) = cellfun(@(x) x-x(1), dataTable(:,3), 'UniformOutput', false);
    idxUnitary = cellfun(@(x) length(x) == 1,dataTable(:,2));
    idxCell    = cellfun('isclass',dataTable(:,2),'cell');
    if (~isempty(find(idxUnitary,1)))
        unitTable = dataTable(idxUnitary & ~idxCell, :);
        unitTable(:,3) = cellfun(@(x) NewTime, unitTable(:,3), 'UniformOutput', false);
        unitTable(:,2) = cellfun(@(x) x*ones(length(NewTime), 1), unitTable(:,2), 'UniformOutput', false);
    end
    
    dataTable = dataTable(~idxUnitary & ~idxCell, :);
%     dataTable(:,2) = cellfun(@(x, y) interp1(double(x), double(y), NewTime, 'linear',y(1)),dataTable(:,3), dataTable(:,2), 'UniformOutput', false);
    dataTable(:,2) = cellfun(@(x,y) interpCellFun(x,y,NewTime),dataTable(:,3), dataTable(:,2), 'UniformOutput', false);
    if (~isempty(find(idxUnitary & ~idxCell,1)))
        dataTable = [dataTable ; unitTable];
    end
    if (SetStartTime0)
        t = NewTime - NewTime(1);
    else
        t = NewTime;
    end
    dataTable = dataTable(:,1:2);
    [~,idx] = sort(upper(dataTable(:,1)));
    finalDatas = [dataTable(idx, :) ; {'t', t}];
end

%% adding for posixtime contained as initialTimestamp of mdfobj
% finalDatas = [finalDatas;{???????'t_posixtime',t+posixtime(mdfObj.InitialTimestamp)}???????];
if Interpolation
    finalDatas = [finalDatas;{'t_posixTime',posixtime(mdfObj.InitialTimestamp)+t}];
end

waitbar(1,curBarH);
close(curBarH)

warning('on','all');

% finalDatas(:,3) = cell(length(finalDatas(:,1)), 1);
if ~isempty(listNames)
    for ii=1:size(finalDatas,1)
        iFirstMatch = 0;
        nameFound = false;
        while ~nameFound && iFirstMatch<size(listNames,1)
            iFirstMatch = iFirstMatch+1;
            nameFound = contains(finalDatas{ii,1},listNames{iFirstMatch});
        end
        if nameFound
            finalDatas{ii,1} = finalDatas{ii,1}(regexp(finalDatas{ii,1},listNames{iFirstMatch},'matchcase','once'):end);
        end
    end
end
finalDatas(:,1) = matlab.lang.makeValidName(finalDatas(:,1));
finalDatas(:,3) = num2cell(ones(length(finalDatas(:,1)) , 1));
end
% FUNCTIONS

% Wildcar function -> filter signals according to listvar
function outputBool = wildcard(sigList,listStatsWith,listEndsWith,listStartAndEndsWith,listMatch)
    startWithBool        = any(cellfun(@(x) any(startsWith(sigList,x)),listStatsWith));
    endWithBool          = any(cellfun(@(x) any(endsWith(sigList,x)),listEndsWith));
    startAndEndsWithBool = any(cellfun(@(x) any(startsWith(sigList,x))& any(endsWith(sigList,x)),listStartAndEndsWith));
    listMatch            = any(cellfun(@(x) any(isequal(sigList,x)),listMatch));
    
    outputBool = startWithBool || endWithBool || startAndEndsWithBool || listMatch;
end
```

Extra-long signals are decoded with the full suffix preserved so that subsequent code can distinguish them:

```
>> who

Your variables are:
            
SUP_DBG_IFDBG_OUT_Output_data_Dev_FUSION3_EgolaneLeftlineLength  
SUP_DBG_IFDBG_OUT_Output_data_Dev_FUSION3_EgolaneLeftline_c0     
SUP_DBG_IFDBG_OUT_Output_data_Dev_FUSION3_EgolaneLeftline_c1     
SUP_DBG_IFDBG_OUT_Output_data_Dev_FUSION3_EgolaneLeftline_psi0   
SUP_DBG_IFDBG_OUT_Output_data_Dev_FUSION3_EgolaneLeftline_x0     
SUP_DBG_IFDBG_OUT_Output_data_Dev_FUSION3_EgolaneLeftline_y0     
SUP_DBG_IFDBG_OUT_Output_data_Dev_FUSION3_EgolaneRightline_c0    
SUP_DBG_IFDBG_OUT_Output_data_Dev_FUSION3_EgolaneRightline_c1    
SUP_DBG_IFDBG_OUT_Output_data_Dev_FUSION3_EgolaneRightline_psi0  
SUP_DBG_IFDBG_OUT_Output_data_Dev_FUSION3_EgolaneRightline_x0    
SUP_DBG_IFDBG_OUT_Output_data_Dev_FUSION3_EgolaneRightline_y0    
SUP_DBG_IFDBG_OUT_Output_data_Dev_FUSION3_previousT0obj_ax       
SUP_DBG_IFDBG_OUT_Output_data_Dev_FUSION3_previousT0obj_ay       
SUP_DBG_IFDBG_OUT_Output_data_Dev_FUSION3_previousT0obj_class    
SUP_DBG_IFDBG_OUT_Output_data_Dev_FUSION3_previousT0obj_vx       
SUP_DBG_IFDBG_OUT_Output_data_Dev_FUSION3_previousT0obj_vy       
SUP_DBG_IFDBG_OUT_Output_data_Dev_FUSION3_previousT0obj_x        
SUP_DBG_IFDBG_OUT_Output_data_Dev_FUSION3_previousT0obj_y        
SUP_DBG_IFDBG_OUT_Output_data_Dev_FUSION3_previousT0obj_yaw         
......
```

## Adapting to the signal format of SWEET 420

After selecting the test path, the tool invokes the scripts `inputFormat.m` to preprocess the data in the .mat files. Due to the different signal names in SWEET 420 and the algorithm added during encoding to shorten signal names, many signals, especially radar signals, cannot be correctly recognized. Consequently, by comparing changes in signal names, we have modified these two scripts to recognize the new signals.

### `ldc_Distance_y`

We had the problem as :

```matlab
Reference to non-existent field 'ldc_Distance_y'.
Error in inputFormat (line 217)
    logConverted.fusion_Distance_y  = logConverted.ldc_Distance_y;
```

`inputFormat.m` Version SWEET 200 :

```matlab
iDistanceLatSig = find(contains(logSignals,'V_m_TargetYdist'),1,'first');
if ~isempty(iDistanceLatSig)
    logConverted.ldc_present = 1;
    logConverted.ldc_Distance_y = log.(logSignals{iDistanceLatSig});
end
```

**`inputFormat.m` Version SWEET 420 :**

```matlab
iDistanceLatSig = find(contains(logSignals,'V_m_Distance_Meas'),1,'first');
if ~isempty(iDistanceLatSig)
    logConverted.ldc_present = 1;
    logConverted.ldc_Distance_y = log.(logSignals{iDistanceLatSig});
end
```

### Radar object signals

We got the warning : 

```
Warning: No radar Object found in current log.
```

#### Filtering radar signals

After clicking the `PostTagging` button, the tool invokes the `inputFormat_v2.m` script, which combines the signals from the .mat files into a data structure in `handles` for GUI plotting. 

It is evident that the camera data of SWEET 420 and SWEET 200 are almost identical. We set breakpoints in the `inputFormat_v2.m` script and then reference the data structures of `camObjects` and `camObjectNames` to generate `radObjects` and `radarObjectNames`.

![image-20240530112320956](C:\Users\p126620\OneDrive - Alliance\Bureau\DEV LOG\assets\image-20240530112320956.png)

![image-20240530112359523](C:\Users\p126620\OneDrive - Alliance\Bureau\DEV LOG\assets\image-20240530112359523.png)

`inputFormat_v2.m` Version SWEET 200 :

```matlab
% Create Radar raw objects
radarObjectSignals = sigNames(contains(sigNames,'Radar_Object') & ~contains(sigNames,'Radar_ObjectAcc') & ~contains(sigNames,'Radar_ObjectAeb') & ~contains(sigNames,'Radar_ObjectInfo') & ~endsWith(sigNames,'_Time'));
radarObjectNames     = unique(cellfun(@(x) x(1:14),radarObjectSignals,'UniformOutput',false));
radObjects = struct();
for i=1:size(radarObjectNames,1)
    currObjectSignals = radarObjectSignals(contains(radarObjectSignals,radarObjectNames{i}));
    for j=1:size(currObjectSignals,1)
        if ~isequal(currObjectSignals{j}(21),'_')
            radObjects(i).(currObjectSignals{j}(21:end)) = log.(currObjectSignals{j});
        else
            radObjects(i).(currObjectSignals{j}(22:end)) = log.(currObjectSignals{j});
        end
    end
end
```

**`inputFormat_v2.m` Version SWEET 420 :**

```matlab
% Create Radar raw objects
radarObjectSignals = sigNames(contains(sigNames, 'IRadar') & ~cellfun(@isempty, regexp(sigNames, '^IRadar.*\d{2}$'))); % the signal begin with 'IRadar' and ends with two digits
radarObjectNames = unique(cellfun(@(x) regexprep(x, '^IRadar.*_(\d{2})$', 'IRadar_Object$1'), radarObjectSignals, 'UniformOutput', false));
numObjects = 64; % SWEET420
radObjects(1:numObjects) = struct();
for i = 1:numel(radarObjectSignals)
    signalName = radarObjectSignals{i};
    objectNumberMatch = regexp(signalName, '_(\d{2})$', 'tokens', 'once');
    objectNumber = str2double(objectNumberMatch{1});
    fieldNameMatch = regexp(signalName, '^IRadar(.*?)_', 'tokens', 'once');
    fieldName = fieldNameMatch{1};
    radObjects(objectNumber + 1).(fieldName) = log.(signalName);
end
```

In SWEET 200 version, most of the radar signals are shaped like `Radar_Object**_XXX`, however in SWEET 420 version, most of the radar signals are shaped like `IRadarXXXObject_**`.

- We revised the filtering rules using the `regexp` function for regular expressions, based on the signal naming conventions in SWEET 420, enabling the program to recognize new radar data.
- And we directly initialize an array of structures of size `numObjects`, reducing runtime memory allocation and management overhead.

Then we have `radObjects` :

![image-20240530115025595](C:\Users\p126620\OneDrive - Alliance\Bureau\DEV LOG\assets\image-20240530115025595.png)

And `radarObjectNames` :

```
>> disp (radarObjectNames)
    {'IRadar_Object00'}
    {'IRadar_Object01'}
    {'IRadar_Object02'}
    {'IRadar_Object03'}
    {'IRadar_Object04'}
    {'IRadar_Object05'}
    
	......
	
    {'IRadar_Object63'}
```

#### Regularization of radar signals

However, the `XXX` part of the new naming convention has also been changed. This resulted in them still not being recognized properly in subsequent `findMPObject.m` scripts. So we modified the field names in `camObjects` in the `inputFormat_v2.m` script to accommodate the algorithm of the postTagging tool.

**In the new `inputFormat_v2.m`:**

```matlab
% Create Radar raw objects
radarObjectSignals = sigNames(contains(sigNames, 'IRadar') & ~cellfun(@isempty, regexp(sigNames, '^IRadar.*\d{2}$'))); % the signal begin with 'IRadar' and ends with two digits
radarObjectNames = unique(cellfun(@(x) regexprep(x, '^IRadar.*_(\d{2})$', 'IRadar_Object$1'), radarObjectSignals, 'UniformOutput', false));
numObjects = 64; % SWEET420
radObjects(1:numObjects) = struct();
for i = 1:numel(radarObjectSignals)
    signalName = radarObjectSignals{i};
    objectNumberMatch = regexp(signalName, '_(\d{2})$', 'tokens', 'once');
    objectNumber = str2double(objectNumberMatch{1});
    fieldNameMatch = regexp(signalName, '^IRadar(.*?)_', 'tokens', 'once');
    fieldName = fieldNameMatch{1};
    radObjects(objectNumber + 1).(fieldName) = log.(signalName);
end

% Copy fields with compatible names
for i = 1:numObjects
    if isfield(radObjects(i), 'LongiDistanceObject')
        radObjects(i).LongitudinalDistanceObject = radObjects(i).LongiDistanceObject;
    end
end
for i = 1:numObjects
    if isfield(radObjects(i), 'AbsLongiVelocityObject')
        radObjects(i).AbsoluteLongiVelocityObject = radObjects(i).AbsLongiVelocityObject;
    end
end
for i = 1:numObjects
    if isfield(radObjects(i), 'AbsLongiAccelObject')
        radObjects(i).AbsoluteLongiAccelerationObject = radObjects(i).AbsLongiAccelObject;
    end
end
```

The names of the joined fields include :

`radObjects.LongitudinalDistanceObject`

`radObjects.AbsoluteLongiVelocityObject`

`radObjects.AbsoluteLongiAccelerationObject`

Since many signal names have changed a lot, there is no means for me to know the relationship between them, and there may also be newly added or removed signals. At present, I modified only these three signal names, so that can already meet the execution of the subsequent program and the basic radar signal identification, if the subsequent found that there are missing radar signals, consider here to retrieve the missing signals.

## Development of pre-processing tool

### Monitoring of context video

In the PostTagging tool, engineers have a requirement to simultaneously monitor videos from the front camara (`_FC.avi`) and context camera (`_RC.avi`), overlaying them in a picture-in-picture (PIP) format.

There are several potential solutions I can think of. For instance, modifying the GUI to create a new video window that can play both videos simultaneously. This approach requires solving numerous technical issues, including :

- Significant difficulty in modifying the GUI.
- Synchronization issues due to differing frame rates between the two videos.
- Modifications needed for video control button functions.
- Performance issues when playing both videos simultaneously, especially when `birdview` is enabled.

Consequently, I opted for a more operable and adaptable solution: writing a preprocessing program that merges the two videos into a new single video beforehand, which is then played normally in the PostTagging tool. This preprocessing program can be accomplished using MATLAB and its Toolboxes and can be integrated with the MF4 decoding tool, or directly invoked by the Perception tool without preprocessing.

Features of the `MergeVideos` function:

- Top comments in the script allow for direct invocation of the function; the if statement on line 78 has been commented out and can be uncommented to use a waitbar for a more efficient way by reduce the refresh rate.
- This code iterates over the subfolders in `directoryPath`, selecting files containing `_FC.avi` and `_RC.avi` for each subfolder, representing front and context videos respectively. The output video file name is based on the front video file name but ends with `_CONTEXT.avi`. Warnings are issued if the output file already exists. If input files are missing, the absence are reported.
- The size and position of the PIP can be freely adjusted as needed.
- For synchronizing the two videos in time, the logic used assumes both videos have the same total duration but different frame rates. The new video's frame rate is set to match the front video, and each frame from the context is evenly distributed over the front video frames. Due to the lower frame rate of the context video, some of its frames are used twice, but they align temporally. Given that the refresh rate of vehicle dashboards (about 1.24 fps) is not high, further optimization is unnecessary.
- The inclusion of a `waitbar` can significantly alleviate anxiety during waiting times. The `waitbar` estimates the remaining processing time based on the time taken for frames already processed.

**Source code of `MergeVideos.m`:**

```matlab
%% The function aimed to merge RC video on FC video and generate a CONTEXT video.
%% The input path should contain capsules as subfolders.

% fprintf('\n-------------VIDEOS WILL BE MERGED IN ALL CAPSULES-----------------\n');
% directoryPath = 'C:\Users\p126620\OneDrive - Alliance\Bureau\MATLAB\TEST MERGE VIDEOS'
% directoryPath = uigetdir;
% MergeVideos1(directoryPath);
% function MergeVideos1(directoryPath)

function MergeVideos(directoryPath)
    subFolders = dir(directoryPath);
    isSub = [subFolders(:).isdir];
    subFolderNames = {subFolders(isSub).name}';
    subFolderNames(ismember(subFolderNames,{'.','..'})) = [];

    h = waitbar(0, 'Initializing...', 'Name', 'Merging Videos Progress', ...
                'Color', [0.97, 0.97, 0.97],'WindowStyle', 'normal');

    for idx = 1:length(subFolderNames)
        subFolderPath = fullfile(directoryPath, subFolderNames{idx});
        files = dir(fullfile(subFolderPath, '*.avi'));
        fileNames = {files.name};
        
        % Input file name
        fcFile = contains(fileNames, '_FC.avi');
        rcFile = contains(fileNames, '_RC.avi');
        
        if any(fcFile) && any(rcFile)
            file1 = fileNames{find(fcFile, 1)};
            file2 = fileNames{find(rcFile, 1)};
            
            video1 = VideoReader(fullfile(subFolderPath, file1));
            video2 = VideoReader(fullfile(subFolderPath, file2));
            
            % Output file name
            outputBaseName = regexprep(file1, 'FC\.avi$', 'CONTEXT.avi');
            outputFileName = fullfile(subFolderPath, outputBaseName);
            
            if exist(outputFileName, 'file')
                warning(['File ', outputFileName, ' already exists.']);
            else
                outputVideo = VideoWriter(outputFileName, 'Motion JPEG AVI');
                outputVideo.FrameRate = video1.FrameRate;
                open(outputVideo);
                
                % Set PIP size and position
                pipWidth = round(video2.Width / pi);
                pipHeight = round(video2.Height / pi);
                pipX = 10; 
                pipY = 10;
                
                totalPipFrames = video2.NumFrames;
                increment = video1.NumFrames / totalPipFrames;
                pipFrameIndex = 1;
                lastPipFrame = readFrame(video2);
                
                tic;
                lastPercentUpdate = 0;
                for k = 1:video1.NumFrames
                    frame1 = readFrame(video1);
                    if k >= round(pipFrameIndex)
                        if hasFrame(video2)
                            lastPipFrame = readFrame(video2);
                            pipFrameIndex = pipFrameIndex + increment;
                        end
                    end
                    
                    resizedPipFrame = imresize(lastPipFrame, [pipHeight pipWidth]);
                    frame1(pipY:(pipY + pipHeight - 1), pipX:(pipX + pipWidth - 1), :) = resizedPipFrame;
                    
                    writeVideo(outputVideo, frame1);
                    
                    currentProgress = round(100 * k / video1.NumFrames);
                    elapsedTime = toc;
                    framesPerSecond = k / elapsedTime;
                    estimatedRemainingTime = (video1.NumFrames - k) / framesPerSecond;
                    
%                     if currentProgress > lastPercentUpdate
                        waitbar(k / video1.NumFrames, h, sprintf('Rendering %s: %d%% - Est. Time: %.2f sec', subFolderNames{idx}, currentProgress, estimatedRemainingTime));
                        lastPercentUpdate = currentProgress;
%                     end
                end
                
                close(outputVideo);
                disp(['COMPLETE: ', subFolderNames{idx}, ' in ', num2str(elapsedTime), ' seconds.']);
            end
        else
            missingFiles = '';
            if ~any(fcFile)
                missingFiles = [missingFiles '_FC.avi '];
            end
            if ~any(rcFile)
                missingFiles = [missingFiles '_RC.avi'];
            end
            disp(['Missing files in ', subFolderNames{idx}, ': ', missingFiles]);
        end
    end
    
    close(h);
end
```

The `waitbar` is shown as:

<img src="C:\Users\p126620\OneDrive - Alliance\Bureau\DEV LOG\assets\image-20240530145252112.png" alt="image-20240530145252112" style="zoom:50%;" />

And we got all the warnings and feedbacks:

```
>> run('C:\Users\p126620\OneDrive - Alliance\Bureau\MATLAB\PreprocessTool\MergeVideos.m')
Missing files in Capsule_Sweet200_1: _FC.avi _RC.avi
Missing files in Capsule_Sweet200_2: _RC.avi
Warning: File C:\Users\p126620\OneDrive - Alliance\Bureau\MATLAB\TEST MERGE
VIDEOS\Capsule_Sweet400\DEV_HJBph2_4066_AD1Evo_VIDEOS_20240409_111831_008_CONTEXT.avi already exists. 
> In MergeVideos>MergeVideos1 (line 40)
In MergeVideos (line 7) 
Missing files in Capsule_Sweet420_1: _RC.avi
COMPLETE: Capsule_Sweet420_2 in 38.4606 seconds.
```

The results in the post tagging tool are shown in the figure, which has no impact on performance.

![image-20240530150310697](C:\Users\p126620\OneDrive - Alliance\Bureau\DEV LOG\assets\image-20240530150310697.png)

### Deleting context videos

We designed a script to delete context videos using the same directory traversal method. Initially we intended to facilitate debugging the `MergeVideos` function. I have retained this script because it may prove useful in the future.

**Source code of `DeleteCONTEXTVideos.m`:**

```matlab
function DeleteCONTEXTVideos(directoryPath)
    dirs = dir(directoryPath);
    isSubfolder = [dirs(:).isdir];
    subfolders = {dirs(isSubfolder).name}';
    subfolders(ismember(subfolders, {'.', '..'})) = []; 
    f = waitbar(0, 'Checking...');
    for k = 1:length(subfolders)
        subfolderName = subfolders{k}; 
        subfolderPath = fullfile(directoryPath, subfolderName);
        videoFiles = dir(fullfile(subfolderPath, '*_CONTEXT.avi'));
        if ~isempty(videoFiles)
            for j = 1:length(videoFiles)
                videoFilePath = fullfile(subfolderPath, videoFiles(j).name);
                delete(videoFilePath);
                fprintf('Deleted: %s/%s\n', subfolderName, videoFiles(j).name);
            end
        else
            fprintf('CONTEXT not found in %s\n', subfolderName);
        end
        waitbar(k / length(subfolders), f, sprintf('Checking: %s', subfolderName));
    end
    close(f);
end
```

### Compiling into an executable program

The `MergeVideos` and `DeleteCONTEXTVideos` functions can be compiled into an executable program, allowing them to be used without the need to launch MATLAB and connect a license. This also makes it easier to distribute to engineers. However, development of the executable program was not pursued further because the Toolbox used for MF4 decoding cannot be compiled.

### Integrating into a single `preprocess` Script

The integration includes:
- Using the same rules for directory selection and folder traversal, supporting one-click processing of multiple capsules.
- Checking data integrity to avoid unnecessary and repetitive tasks, providing detailed feedback.
- Addressing issues with the progress bar and processing order.
- Solve the problem of repeated calls to external functions and data, the script can now run independently without relying on other files.

**Source code of `preprocess.m`:**

> [!NOTE]
>
> The code is too long so we put it at the end of the document.

### User experience
Despite not being able to compile to an executable program, I made a simple command line UI that is easier to use and reduces the learning cost of the tool for engineers.

```matlab
    while true
        clc;
        disp('====================================');
        disp(' Preprocess Tool for Perception ');
        disp('====================================');
        disp('[1] Preprocess MF4 and Videos');
        disp('[2] Only Generate Videos');
        disp('[3] Delete Videos');
        disp('[4] Exit');
        disp('====================================');
        choice = input('Enter your choice : ');

        switch choice
            case 1
                clc;
                disp('====================================');
                disp(' Preprocess MF4 and Videos ');
                disp('====================================');
                run('preprocess.m')
                pause(0.5);
                    confirm2 = input('Continue ? (y/n): ', 's');
                    if lower(confirm2) == 'y'
                    else
                    end
            case 2
                directoryPath = uigetdir;
                fprintf('\n--------VIDEOS WILL BE MERGED IN ALL CAPSULES SELECTED------------\n');
                MergeVideos(directoryPath);
                pause(5);
            case 3
                directoryPath = uigetdir;
                fprintf('\n--------ALL THE MERGED VIDEOS SELECTED WILL BE REMOVED------------\n');
                DeleteCONTEXTVideos(directoryPath);
                pause(5);
            case 4
                clear all;
                break;
            otherwise
        end
    end
```

When we run the tool:

```
====================================
 Preprocess Tool for Perception 
====================================
[1] Preprocess MF4 and Videos
[2] Only Generate Videos
[3] Delete Videos
[4] Exit
====================================
Enter your choice : 
```

We select the first choice and we select a path which has four capsules for test:

```
====================================
 Preprocess MF4 and Videos 
====================================
 Capsules found : 
    {'Capsule_Sweet200_1'}
    {'Capsule_Sweet400'  }
    {'Capsule_Sweet420_1'}
    {'Capsule_Sweet420_2'}

1/4 : Creating 20210422_103639.mat log.
	1/5 : Converting _CAN-ITS1_ : DEV_HHN_626_AD1Evo_CAN-ITS1_20210422_103641_001_ESW-FCam7.2-FRad7.1.2.MF4	 -> Done !
	2/5 : Converting _CAN-ITS2_ : DEV_HHN_626_AD1Evo_CAN-ITS2_20210422_103641_001_ESW-FCam7.2-FRad7.1.2.MF4	 -> Done !
	3/5 : Converting _CAN-ITS4_ : DEV_HHN_626_AD1Evo_CAN-ITS4_20210422_103641_001_ESW-FCam7.2-FRad7.1.2.MF4	 -> Done !
	4/5 : Converting _FCam_ : DEV_HHN_626_AD1Evo_FCam_20210422_103641_001_ESW-FCam7.2-FRad7.1.2.MF4	 -> Done !
	5/5 : Converting _FRad_ : DEV_HHN_626_AD1Evo_FRad_20210422_103641_001_ESW-FCam7.2-FRad7.1.2.MF4	 -> Done !
		  Converting VIDEO : DEV_HHN_626_AD1Evo_VIDEOS_20210422_103641_001.MF4	 -> Done !
	      Saving 20210422_103641.mat...	 -> Done !
          Missing files in Capsule_Sweet200_1: _FC.avi _RC.avi

2/4 : Creating 20240409_111828.mat log.
	1/5 : Converting _CAN-ITS1_ : DEV_HJBph2_4066_AD1Evo_CAN-ITS1_20240409_111830_008_ESW-FCam3.0-FRad3.3.MF4	 -> Done !
	2/5 : Converting _CAN-ITS2_ : DEV_HJBph2_4066_AD1Evo_CAN-ITS2_20240409_111830_008_ESW-FCam3.0-FRad3.3.MF4	 -> Done !
	3/5 : Converting _CAN-ITS4_ : DEV_HJBph2_4066_AD1Evo_CAN-ITS4_20240409_111830_008_ESW-FCam3.0-FRad3.3.MF4	 -> Done !
	4/5 : Converting _FCam_ : DEV_HJBph2_4066_AD1Evo_FCam_20240409_111830_008_ESW-FCam3.0-FRad3.3.MF4	 -> Done !
	5/5 : Converting _FRad_ : DEV_HJBph2_4066_AD1Evo_FRad_20240409_111830_008_ESW-FCam3.0-FRad3.3.MF4	 -> Done !
		  Converting VIDEO : DEV_HJBph2_4066_AD1Evo_VIDEOS_20240409_111831_008.MF4	 -> Done !
	      Saving 20240409_111830.mat...	 -> Done !
          VIDEO MERGING COMPLETE: Capsule_Sweet400

3/4 : Creating 20240430_092642.mat log.
	1/5 : Converting _CAN-ITS1_ : DEV_HJBph2_4066_AD1Evo_CAN-ITS1_20240430_092644_007_ESW-FCam3.0-FRad3.3.MF4	 -> Done !
	2/5 : Converting _CAN-ITS2_ : DEV_HJBph2_4066_AD1Evo_CAN-ITS2_20240430_092644_007_ESW-FCam3.0-FRad3.3.MF4	 -> Done !
	3/5 : Converting _CAN-ITS4_ : DEV_HJBph2_4066_AD1Evo_CAN-ITS4_20240430_092644_007_ESW-FCam3.0-FRad3.3.MF4	 -> Done !
	4/5 : Converting _FCam_ : DEV_HJBph2_4066_AD1Evo_FCam_20240430_092644_007_ESW-FCam3.0-FRad3.3.MF4	 -> Done !
          !! No recorder _FRad_ found at 20240430_092642 !! 
		  Converting VIDEO : DEV_HJBph2_4066_AD1Evo_VIDEOS_20240430_092645_007.MF4	 -> Done !
	      Saving 20240430_092644.mat...	 -> Done !
          Missing files in Capsule_Sweet420_1: _RC.avi

4/4 : Creating 20240430_140635.mat log.
	1/5 : Converting _CAN-ITS1_ : DEV_HJBph2_4066_AD1Evo_CAN-ITS1_20240430_140637_007_ESW-FCam3.0-FRad3.3.MF4	 -> Done !
	2/5 : Converting _CAN-ITS2_ : DEV_HJBph2_4066_AD1Evo_CAN-ITS2_20240430_140637_007_ESW-FCam3.0-FRad3.3.MF4	 -> Done !
	3/5 : Converting _CAN-ITS4_ : DEV_HJBph2_4066_AD1Evo_CAN-ITS4_20240430_140637_007_ESW-FCam3.0-FRad3.3.MF4	 -> Done !
	4/5 : Converting _FCam_ : DEV_HJBph2_4066_AD1Evo_FCam_20240430_140637_007_ESW-FCam3.0-FRad3.3.MF4	 -> Done !
	5/5 : Converting _FRad_ : DEV_HJBph2_4066_AD1Evo_FRad_20240430_140637_007_ESW-FCam3.0-FRad3.3.MF4	 -> Done !
		  Converting VIDEO : DEV_HJBph2_4066_AD1Evo_VIDEOS_20240430_140638_007.MF4	 -> Done !
	      Saving 20240430_140637.mat...	 -> Done !
          VIDEO MERGING COMPLETE: Capsule_Sweet420_2

Continue ? (y/n): 
```

Thus we can see all the detailed logs before continue using the post tagging tool.

##  Plotting multiple signals

In the PostTagging tool, engineers wish to simultaneously plot multiple signals within the drawing area and view legends. This task presents significant challenges, and I aim to implement this feature in a straightforward manner, starting by modifying the `refreshGraph.m` script. This script is continually invoked when the user interacts with the PostTagging GUI, used to update data in `handles` and redraw the plot. My approach involves inserting a `hold on` command at the appropriate position within the script and setting up a `legend`.

```matlab
function handles = refreshGraph(handles, onlyMarker)
    if nargin < 2 
        onlyMarker = 0;
    end
    global currSignal
    global currTime
    plotValue = handles.loadedLog.(currSignal);
    
    if size(handles.loadedLog.t, 1) == size(plotValue, 1)
        currentLogTime = round(currTime-handles.vCont.t0+handles.loadedLog.t(1), 3);
        [tUnique, iTimeUnique, ~] = unique(handles.loadedLog.t);
        currIndexe = ceil(interp1(tUnique, iTimeUnique, currentLogTime));
        
        axes(handles.Graph); 
        if ~isnan(currIndexe)
            if onlyMarker == 0
                randomColor = rand(1,3); 
                hold on; 
                    hLine = plot(handles.loadedLog.t, plotValue, 'Color', randomColor, 'DisplayName', sprintf('%s Signal', currSignal));
                hold off;
            end
            set(handles.markerPlot, 'XData', handles.loadedLog.t(currIndexe), 'YData', plotValue(currIndexe), 'DisplayName', 'Current Marker');
        end
        
        if isequal(currSignal, 'TAG_ACC')
            if ~any(contains(handles.Graph.YTickLabel, 'RAS'))
                set(handles.Graph, 'YTick', [0 1 2 3 4 5 6 10], 'YTickLabel', {'RAS', 'LS', 'LD', 'NLD', 'TOL', 'GHOST', 'ACC', 'Other'});
                
                iTagACC_K1_val = find(round(rem(plotValue, 1), 1) == 0.1 & plotValue ~= 0);
                iTagACC_K2_val = find(round(rem(plotValue, 1), 1) == 0.2 & plotValue ~= 0);
                iTagACC_K3_val = find(round(rem(plotValue, 1), 1) == 0.3 & plotValue ~= 0);
                
                
                hold on; 
                plot(handles.loadedLog.t(iTagACC_K1_val), plotValue(iTagACC_K1_val), 'g*', 'DisplayName', 'TAG_ACC K1');
                plot(handles.loadedLog.t(iTagACC_K2_val), plotValue(iTagACC_K2_val), 'm*', 'DisplayName', 'TAG_ACC K2');
                plot(handles.loadedLog.t(iTagACC_K3_val), plotValue(iTagACC_K3_val), 'y*', 'DisplayName', 'TAG_ACC K3');
                hold off; 
            end
        else
            if any(contains(handles.Graph.YTickLabel, 'RAS'))
                handles.Graph.YTickLabelMode = 'Auto';
                handles.Graph.YTickMode = 'Auto';
            end
        end
        legend(handles.Graph, 'AutoUpdate', 'on', 'TextColor', 'white'); 
    else
        error('Wrong signal selected!');
    end
end
```

**Source code of `preprocess.m`:**

```matlab
% clear all
% fprintf('\n-------------Preprocess MF4 and Videos-----------------\n');
% directoryPath = 'C:\Users\p126620\OneDrive - Alliance\Bureau\MATLAB\TEST MERGE VIDEOS'
    directoryPath = uigetdir; % select the path with subfolders which contain MF4 and AVI files
%   toolbox functions required
 	preprocess(directoryPath);
    subFolders = dir(directoryPath);
    isSub = [subFolders(:).isdir];
    subFolderNames = {subFolders(isSub).name}';
    subFolderNames(ismember(subFolderNames,{'.','..'})) = [];
    disp(' Capsules found : ');
    disp(subFolderNames);

    h = waitbar(0, 'Initializing...', 'Name', 'Merging Videos Progress', ...
                'Color', [0.97, 0.97, 0.97],'WindowStyle', 'normal');

    for idx = 1:length(subFolderNames)
        subFolderPath = fullfile(directoryPath, subFolderNames{idx});
        
        %% ConvertLogs
        maxTime             = 90;
        signalBased         = 1; %1 : Convert only signal based records
                                 %0 : Convert only message based records
        concatenateLogsAnswer = 0;
        exportTXTFile       = 0;
        % {recorderName recordingType (signalBased=1/messageBased=0)}
        recorders = {
        %              '_CAN-ADAS_'
        %              '_Adas_'
                     '_CAN-ITS1_'   1
                     '_CAN-ITS2_'   1
        %              '_CAN-ITS3_'   0
                     '_CAN-ITS4_'   1
                     '_FCam_'       1
                     '_FRad_'       1
                    };
        count_its2 = 0;
        iDataFrame= zeros(size(recorders,1),1);

                                    %% CANAPE LOGS

                                    listVarCAN = {
                                                  % ADAS ecu
                                    %               '*_ITS2'
                                  'V_*'
                                  'TP_*'
                                  'F_x*'
                                  'ISteeringWheelAngleCorrected*'

                                  % GT
                                  'device_1_Analog_0' % GT button


                                  % Target flag
                                  '*RadarIDobject' % Bosch Radar
                                  '*RadarCIPVobject'
                                  '*CamCIPVobject' % CIPV Camera

                                  % Longitudinal Distance
                                 '*LongitudinalDistanceObject'
                                 'Range1PosForward' % RT Range

                                 % Lateral Distance
                                 '*LateralDistanceObject'
                                 'Range1PosLateral' 

                                 % Lateral Speed
                                 '*RelativeLatVelocityObject'
                                 'Range1VelLateral'

                                 % Lateral Accel
                                 '*RelativeLatAccelerationObject'

                                 % Longitudinal Speed
                                 '*AbsoluteLongiVelocityObject'
                                 'Range1VelForward'
                                 'HunterVelForward'

                                 % Longitudinal Accel
                                 '*AbsoluteLongiAccelerationObject'

                                 % Time stamps
                                 '*RadarTimestampObject' % Bosch Radar

                                 % Other
                                 'externalClockTimestamp'
                                 'TAG_ACC'

                                 % Other object infos
                                 '*HeadingAngleObject'
                                 '*CamAzimuthObject'
                                 '*LengthObject'
                                 '*WidthObject'
                                 '*IDobject'
                                 '*ClassificationObject'
                                 '*AzimuthObject'
                                 '*ExistenceProbObject'

                                 % Lines
                                 % Left Left
                                 'Cam_InfrastructureLines_CamLeftLeftLine*'
                                 % Left
                                 'Cam_InfrastructureLines_CamLeftLine*'
                                 % Right
                                 'Cam_InfrastructureLines_CamRightLine*'
                                 % Right Right
                                 'Cam_InfrastructureLines_CamRightRightLine*'

                                 % Other Range
                                 'Range1TimeToCollisionForward'
                                 'Range1TimeToCollisionLateral'
                                 'Range1VelLateral'

                                 };


                                    %% CONTROL DESK LOGS
                                    listVarCD = {% HHN FUsion Object 1
                                                'TP_HHN_1E'
                                                'TP_HHN_2E'

                                                % Longitudinal Distance
                                                'TP_V_m_Distance' % Fusion RSA
                                                'TP_C1A_Vector_RadarXxObject' % Radar Conti
                                                'TP_C1A_Vector_RadarXxObject01' % Radar Conti
                                                'TP_C1A_Vector_RadarXxObject00XxDistanceObject00' % Radar Conti
                                                'TP_C1A_Vector_CameraXxCVM_Object' % Camera Valeo
                                                'TP_Range1PosForward' % RT Range

                                                % Taget Flag
                                                'F_x_Target' % Fusion RSA
                                                'TP_C1A_Vector_RadarXxObject00XxIDobject00' % Radar Conti
                                                'TP_C1A_Vector_CameraXxCVM_Object00XxCamCIPVobject00' % Camera Valeo

                                                % Lateral Distance
                                                'V_m_TargetYdistCalc' % Fusion RSA
                                                'TP_C1A_Vector_RadarXxObject00XxLateralDistanceObject00' % Radar Conti
                                                'TP_C1A_Vector_CameraXxCVM_Object00XxCamLateralDistanceObject00' % Camera Valeo
                                                'TP_Range1PosLateral' % RT Range

                                                % Longitudinal Speed
                                                'V_mps_TargetSpeed' % Fusion RSA
                                                'TP_C1A_Vector_RadarXxObject00XxRelativeSpeedObject00' % Radar Conti
                                                'TP_C1A_Vector_CameraXxCVM_Object00XxCamRelativeVelocityObject00' % Camera Valeo
                                                'TP_Range1VelForward' % RT Range

                                                % Longitudinal Accel
                                                'TP_V_mps2_TargetAcceleration_Est' % Fusion RSA
                                                'TP_C1A_Vector_RadarXxObject00XxRelativeAccelerationObject00' % Radar Conti
                                                'TP_C1A_Vector_CameraXxCVM_Object00XxCamAbsoluteAccelerationObject00' % Camera Valeo

                                                % Other
                                                'externalclocktimestamp' % external clock time stamp for synchro
                                                'TP_V_mes_mps_VehicleSpeed' % Ego Speed
                                                'TP_V_mps2_VehAcceleration_Est' % Ego Acceleration
                                                'V_degps_YawRate' % Ego Yaw Rate
                                                'V_deg_SteeringWheelAngle' % Steering Wheel Angle
                                                'TP_F_x_'
                                                '=t'

                                                % LDC Fusion selection
                                                'TP_Flag_Fusion_P_tng_x_AD1Enh_FusionFeatureMngt' % C1AHS Fusion used Flag
                                                }';

                                       % Old
                                    %    listVarCAN = {
                                    %               % ADAS ecu
                                    %               'Target'
                                    %               'V_mes_mps_VehicleSpeed_ITS2'
                                    %               'V_mps2_VehicleAcceleration_Est_ITS2'
                                    %               'IYawRateCorrected'
                                    %               'ISteeringWheelAngleCorrected'
                                    %               
                                    %               % ADAS dev frames
                                    %               'F_x_TargetDetected_ITS2'
                                    %               'V_m_Distance_Meas_ITS2'
                                    %               'V_mps_TargetSpeed_Meas_ITS2'
                                    %               'V_mps2_TargetAcceleration_ITS2'
                                    %               'V_m_TargetYdistMeas_ITS2'
                                    %               
                                    %               
                                    %               
                                    %               'TP_V_mes_mps_VehicleSpeed' % Ego Speed
                                    %               'TP_V_mps2_VehAcceleration_Est' % Ego Acceleration
                                    %               'V_degps_YawRate' % Ego Yaw Rate
                                    %               'V_deg_SteeringWheelAngle' % Steering Wheel Angle
                                    %               
                                    %             
                                    %               % GT
                                    %               'device_1_Analog_0' % GT button
                                    %               
                                    % 
                                    %               % Target flag
                                    %               'Radar_ObjectAcc00_RadarIDobject' % Bosch Radar
                                    %               'Cam_Object00_CamCIPVobject' % ZF Camera
                                    %               'Radar_ObjectAcc01_RadarIDobject' % Bosch Radar
                                    %               'Cam_Object01_CamCIPVobject' % ZF Camera
                                    %                 
                                    %               % Longitudinal Distance
                                    %              'Radar_ObjectAcc00_RadarLongitudinalDistanceObject' % Bosch Radar
                                    %              'Cam_Object00_CamLongitudinalDistanceObject' % ZF Camera
                                    %              'Radar_ObjectAcc01_RadarLongitudinalDistanceObject' % Bosch Radar
                                    %              'Cam_Object01_CamLongitudinalDistanceObject' % ZF Camera
                                    %              'Range1PosForward' % RT Range
                                    %              
                                    %              % Lateral Distance
                                    %              'Radar_ObjectAcc00_RadarLateralDistanceObject' % Bosch Radar
                                    %              'Cam_Object00_CamLateralDistanceObject' % ZF Camera
                                    %              'Radar_ObjectAcc01_RadarLateralDistanceObject' % Bosch Radar
                                    %              'Cam_Object01_CamLateralDistanceObject' % ZF Camera
                                    %              'Range1PosLateral' 
                                    %              
                                    %              % Longitudinal Speed
                                    %              'Radar_ObjectAcc00_RadarAbsoluteLongiVelocityObject' % Bosch Radar
                                    %              'Cam_Object00_CamAbsoluteLongiVelocityObject' % ZF Camera
                                    %              'Radar_ObjectAcc01_RadarAbsoluteLongiVelocityObject' % Bosch Radar
                                    %              'Cam_Object01_CamAbsoluteLongiVelocityObject' % ZF Camera
                                    %              'Range1VelForward'
                                    %              'HunterVelForward'
                                    %              
                                    %              % Longitudinal Accel
                                    %              'Radar_ObjectAcc00_RadarAbsoluteLongiAccelerationObject' % Bosch Radar
                                    %              'Cam_Object00_CamAbsoluteLongiAccelerationObject' % ZF Camera
                                    %              'Radar_ObjectAcc01_RadarAbsoluteLongiAccelerationObject' % Bosch Radar
                                    %              'Cam_Object01_CamAbsoluteLongiAccelerationObject' % ZF Camera
                                    %              
                                    %              % Time stamps
                                    %              'Radar_ObjectAcc00_RadarTimestampObject' % Bosch Radar
                                    %              'Cam_Object00_CamTimestampObject' % ZF Camera
                                    %              'Radar_ObjectAcc01_RadarTimestampObject' % Bosch Radar
                                    %              'Cam_Object01_CamTimestampObject' % ZF Camera
                                    %              
                                    %              % Other
                                    %              'externalClockTimestamp'
                                    %              'TAG_ACC'
                                    %              
                                    %              % Other object infos
                                    %              'Cam_Object00_CamHeadingAngleObject'
                                    %              'Cam_Object01_CamHeadingAngleObject'
                                    %              'Cam_Object00_CamAzimuthObject'
                                    %              'Cam_Object01_CamAzimuthObject'
                                    %              'Cam_Object00_CamLengthObject'
                                    %              'Cam_Object01_CamLengthObject'
                                    %              'Cam_Object00_CamWidthObject'
                                    %              'Cam_Object01_CamWidthObject'
                                    %              
                                    %              % Lines
                                    %              % Left Left
                                    %              'Cam_InfrastructureLines_CamLeftLeftLineOffset'
                                    %              'Cam_InfrastructureLines_CamLeftLeftLineYawAngle'
                                    %              'Cam_InfrastructureLines_CamLeftLeftLineCurvature'
                                    %              'Cam_InfrastructureLines_CamLeftLeftLineCurvatureRate'
                                    %              % Left
                                    %              'Cam_InfrastructureLines_CamLeftLineOffset'
                                    %              'Cam_InfrastructureLines_CamLeftLineYawAngle'
                                    %              'Cam_InfrastructureLines_CamLeftLineCurvature'
                                    %              'Cam_InfrastructureLines_CamLeftLineCurvatureRate'
                                    %              % Right
                                    %              'Cam_InfrastructureLines_CamRightLineOffset'
                                    %              'Cam_InfrastructureLines_CamRightLineYawAngle'
                                    %              'Cam_InfrastructureLines_CamRightLineCurvature'
                                    %              'Cam_InfrastructureLines_CamRightLineCurvatureRate'
                                    %              % Right Right
                                    %              'Cam_InfrastructureLines_CamRightLineOffset'
                                    %              'Cam_InfrastructureLines_CamRightLineYawAngle'
                                    %              'Cam_InfrastructureLines_CamRightLineCurvature'
                                    %              'Cam_InfrastructureLines_CamRightLineCurvatureRate'
                                    %              
                                    %              % Other Range
                                    %              'Range1TimeToCollisionForward'
                                    %              'Range1TimeToCollisionLateral'
                                    %              'Range1VelLateral'
                                    %              
                                    %              % Driver Brake
                                    %              'IBrakePedalPressedByDriver'
                                    %              'F_x_DriverBrake_ITS2'
                                    %              };
        listVarCAN = {};
        currCanapePath = subFolderPath;
            % canapeConvPath = fullfile(currCanapePath,'test');
            % Convert CANape Log
            canapeFiles = filesearch(currCanapePath,'MF4',0);
            canapeNamesRaw = {canapeFiles.name}';
            canapeNames = canapeNamesRaw(~contains(canapeNamesRaw,'zf_frCam.MF4') & ~contains(canapeNamesRaw,'RadarFC.MF4') & ~contains(canapeNamesRaw,'RIF_FC.MF4'));

            videoFiles  = filesearch(currCanapePath,'AVI',0);
            videoNames  = {videoFiles.name}';
            videoFCNames= videoNames(contains(videoNames,'_FC.avi'));

            matFiles                = dir([currCanapePath '\*.mat']);
            matNamesAll             = {matFiles.name}';
            matNamesCapsules        = {};
            matNamesConcatenations  = {};
            if ~isempty(matNamesAll)
                iCapsules           = find(~contains(matNamesAll,'_total.mat'));
                iConcatenations     = find(contains(matNamesAll,'_total.mat'));
                if ~isempty(iCapsules)
                    matNamesCapsules    = matNamesAll(iCapsules);
                end
                if ~isempty(iConcatenations)
                    matNamesConcatenations = matNamesAll(iConcatenations);
                end
            end


            fileNames = canapeNames(contains(canapeNames,recorders(:,1)));
            videosFiles= canapeNames(contains(canapeNames,'VIDEOS'));
            gpsFiles    = canapeNames(contains(canapeNames,'_GPSLabelling_'));

            timeCharArray = {};
            timePosixArray = [];
            for f=1:size(fileNames,1)
                timeCharacter = regexp(fileNames{f},'\d{8}_\d{6}','match','once');
                Y   = str2num(timeCharacter(1:4));
                M   = str2num(timeCharacter(5:6));
                D   = str2num(timeCharacter(7:8));
                H   = str2num(timeCharacter(10:11));
                Mi  = str2num(timeCharacter(12:13));
                S   = str2num(timeCharacter(14:15));
                timeDate = datetime(Y,M,D,H,Mi,S);
                timePosix    = posixtime(timeDate);
                iCloseFile   = find(abs(timePosixArray-timePosix)<=2); % find all previous record 2s neir current one
                if isempty(iCloseFile) % if no record close to current one already stored
                    potentialDateTimes  = datetime([timePosix-2;timePosix-1;timePosix;timePosix+1;timePosix+2], 'ConvertFrom', 'posixtime');
                    potentialTimeCharacters = cellstr(datestr(potentialDateTimes,'yyyymmdd_HHMMSS'))';
                    timeCharArray = [timeCharArray;potentialTimeCharacters];
                    timePosixArray = [timePosixArray;timePosix];
                end
            end
        %     tic
            for i = 1:size(timeCharArray,1)
                if any(ismember(matNamesCapsules,[timeCharArray{i,3} '.mat'])) % mat file converted already exists
        %             eraseMatFileAnswer = questdlg(sprintf('%s already exists. Ovewrite it ?',[timeCharArray{i,3} '.mat']),'.mat already exists','Ignore','Overwrite','Ignore');
                    if false %isequal(eraseMatFileAnswer,'Ignore')
                        fprintf('\n%d/%d : %s.mat already exists, Ignored.\n',i,size(timeCharArray,1),timeCharArray{i});
                        continue;
                    end
                end
                % Check if mat exists
                outputFileName = fullfile(currCanapePath, [timeCharArray{i,3} '.mat']);
                if exist(outputFileName, 'file')
                    fprintf('%s already exists. Skipping processing.\n', outputFileName);
                else
                
                capsuleFiles = fileNames(contains(fileNames,timeCharArray(i,:)));
                fprintf('\n%d/%d : Creating %s.mat log.\n',i,size(timeCharArray,1),timeCharArray{i});
                clear logMergedCropped logVideo logVideoCropped
                logMerged = struct();
                count_its2 = 0;
                for rec = 1:size(recorders,1)
                    clear log logCropped logMergedCropped logVideo log

                    iCorrespLog = find(contains(capsuleFiles,recorders{rec,1}));
                    if ~isempty(iCorrespLog) % if corresponding .MF4 found
                        recorderFile = capsuleFiles{iCorrespLog};
                        if recorders{rec,2} == 1 % Signal based conversion.
                            fprintf('\t%d/%d : Converting %s : %s',rec,size(recorders,1),recorders{rec,1},recorderFile);
                            log = readMDF_DM_MD(fullfile(currCanapePath,recorderFile), 1, 0.01, 0,listVarCAN);
                            fprintf('\t -> Done !\n');
                        else % messageBased conversion
                           if contains(recorders{rec,1},'FCam') % FCam recorder
                               fprintf('\t%d/%d : Converting FCam : %s',rec,size(recorders,1),recorderFile);
                               [log iDataFrame(rec)]   = readMDF_MsgBased_Ethernet(fullfile(currCanapePath,recorderFile), 1, 0.01, 1,1,iDataFrame(rec),listVarCAN);
                               fprintf('\t -> Done !\n');
            %                    toc
                           elseif contains(recorders{rec,1},'FRad') % FRad recorder
                               fprintf('\t%d/%d : Converting FRad : %s',rec,size(recorders,1),recorderFile);
                               [log iDataFrame(rec)]   = readMDF_MsgBased_Ethernet(fullfile(currCanapePath,recorderFile), 1, 0.01, 1,2,iDataFrame(rec),listVarCAN);
                               log = compuCoefs(log,'RadarLongitudinalDistanceObject',0.078125,1);
                               log = compuCoefs(log,'RadarLateralDistanceObject',0.078125,1);
                               log = compuCoefs(log,'RadarAbsoluteLongiVelocityObject',0.00390625,1);
                               log = compuCoefs(log,'RadarAbsoluteLatVelocityObject',0.00390625,1);
                               log = compuCoefs(log,'RadarRelativeLatVelocityObject',0.00390625,1);
                               log = compuCoefs(log,'RadarAbsoluteLongiAccelerationObject',0.00048828125,1);
                               log = compuCoefs(log,'RadarAbsoluteLatAccelerationObject',0.00048828125,1);
                               log = compuCoefs(log,'RadarRelativeLatAccelerationObject',0.00048828125,1);
                               log = compuCoefs(log,'RadarLengthObject',0.0078125,1);
                               log = compuCoefs(log,'RadarWidthObject',0.0078125,1);
                               log = compuCoefs(log,'RadarHeadingAngleObject',-0.99988,1);
                               fprintf('\t -> Done !\n');
                           elseif contains(recorders{rec,1},'CAN-ADAS') || (contains(recorders{rec,1},'_CAN-ITS2_') && count_its2 == 1) % CAN-ADAS -> TVC recorder.
                               fprintf('\t%d/%d : Converting ADAS : %s',rec,size(recorders,1),recorderFile);
                               log          = readMDF_DM_MD(fullfile(currCanapePath,recorderFile), 1, 0.01, 0); % 
                               if isstruct(log) && ~isempty(fieldnames(log))
                                    fprintf('\t -> Done !\n');
                               end
                           elseif contains(recorders{rec,1},'Adas') % Adas -> ADAS internal signals.
                               fprintf('\t%d/%d : Converting Adas : %s',rec,size(recorders,1),recorderFile);
                               log          = readMDF_DM_MD(fullfile(currCanapePath,recorderFile), 1, 0.01, 0,{},{'adegocenterlane','adegooal','adlineoal','adtargetoal'}'); % 

                               if isstruct(log) && ~isempty(fieldnames(log))
                                    fprintf('\t -> Done !\n');
                               end
            %                    toc
                           else% its
                               if contains(recorders{rec,1},'ITS1')
                                   fprintf('\t%d/%d : Converting ITS1 : %s',rec,size(recorders,1),recorderFile);
                                   itsNb = 1;
                               elseif contains(recorders{rec,1},'ITS2')
                                   fprintf('\t%d/%d : Converting ITS2 : %s',rec,size(recorders,1),recorderFile);
                                   itsNb = 2;
                                   count_its2 = 1;
                               elseif contains(recorders{rec,1},'ITS3')
                                   fprintf('\t%d/%d : Converting ITS3 : %s',rec,size(recorders,1),recorderFile);
                                   itsNb = 3;
                               elseif contains(recorders{rec,1},'ITS4')
                                   fprintf('\t%d/%d : Converting ITS4 : %s',rec,size(recorders,1),recorderFile);
                                   itsNb = 4;
                               elseif contains(recorders{rec,1},'ITS5')
                                   fprintf('\t%d/%d : Converting ITS5 : %s',rec,size(recorders,1),recorderFile);
                                   itsNb = 5;
                               end
                               [log iDataFrame(rec)]   = readMDF_MsgBased_CAN_FD(fullfile(currCanapePath,recorderFile), 1, 0.01, 1,itsNb,iDataFrame(rec),listVarCAN);
                               fprintf('\t -> Done !\n');
            %                    toc
                           end
                        end
                        if ~isempty(fieldnames(log)) && ~isempty(log.t)
                            log.t           = log.t - log.t(1);
                            if floor(log.t(end)) > 90
                                iOk             = find(log.t > floor((log.t(end)-0.1)/90)*90);
                                log             = cropLog(log,iOk(1),iOk(end),100);
                            end


                            if isempty(fieldnames(logMerged))
                                logMerged = log;
                            else
                                minPosixTime = max(log.t_posixTime(1),logMerged.t_posixTime(1));
                                maxPosixTime = min(log.t_posixTime(end),logMerged.t_posixTime(end));

                                [closestSample begin1] = min(abs(log.t_posixTime-minPosixTime));
                                [closestSample end1]   = min(abs(log.t_posixTime-maxPosixTime));
                                [closestSample begin2] = min(abs(logMerged.t_posixTime-minPosixTime));
                                [closestSample end2]   = min(abs(logMerged.t_posixTime-maxPosixTime));

                                if end1-begin1 ~= end2-begin2
                                    if end1-begin1 == end2-begin2+1
                                        end1 = end1-1;
                                    elseif end1-begin1 == end2-begin2-1
                                        begin1   = begin1-1;
                                    else
                                        error('Logs have more than sample time delay !');
                                    end
                                end 
            %                     [correlFound,begin1,end1,begin2,end2] = findCorrelation(log.IYawRateCorrected,logMerged.IYawRateCorrected,100,100);
                                if end1-begin1 >0 % at least 1 sample time in common
                                    logCropped = cropLog(log,max(begin1,1),min(size(log.t,1),end1),100);
                                    if begin2>1 || end2<size(logMerged.t,1)
                                        logMerged = cropLog(logMerged,max(begin2,1),min(size(logMerged.t,1),end2),100);
                                    end
                                    logMerged = appendStructs(logMerged,logCropped);
                                else
                                    fprintf('\n No correlation found at %s with file %s !! \n',timeCharArray{i},recorderFile);
                                end
                            end
                        end
                    else
                        fprintf('\n !! No recorder %s found at %s !! \n',recorders{rec,1},timeCharArray{i});
                    end
                end

                if ~isempty(fieldnames(logMerged)) & ~isempty(logMerged.t)
                    % Load video Recorder
                    if ~isempty(videosFiles) && any(contains(videosFiles,timeCharArray(i,:)))
                        recorderVideoFile = videosFiles{contains(videosFiles,timeCharArray(i,:))};
                        fprintf('\t\t Converting VIDEO : %s',recorderVideoFile);
                        logVideo      = readMDF_DM_MD(fullfile(currCanapePath,recorderVideoFile), 1, 0.01, 0);
                        if ~isempty(fieldnames(logVideo))
                            fprintf('\t -> Done !\n');
                            if ~isempty(logVideo.t)
                %             toc
                                minPosixTime = max(logVideo.t_posixTime(1),logMerged.t_posixTime(1));
                                maxPosixTime = min(logVideo.t_posixTime(end),logMerged.t_posixTime(end));

                                [closestSample iBeginVideo] = min(abs(logVideo.t_posixTime-minPosixTime));
                                [closestSample iEndVideo]   = min(abs(logVideo.t_posixTime-maxPosixTime));
                                [closestSample iBeginLog] = min(abs(logMerged.t_posixTime-minPosixTime));
                                [closestSample iEndLog]   = min(abs(logMerged.t_posixTime-maxPosixTime));
                                if iEndVideo-iBeginVideo ~= iEndLog-iBeginLog
                                    if iEndVideo-iBeginVideo == iEndLog-iBeginLog+1
                                        iEndVideo = iEndVideo-1;
                                    elseif iEndVideo-iBeginVideo == iEndLog-iBeginLog-1
                                        iBeginVideo = iBeginVideo-1;
                                    else
                                        error('Video and Log have more than sample time delay !');
                                    end
                                end
                                if iBeginLog>1 || iEndLog<size(logMerged.t,1)
                                    logMerged = cropLog(logMerged,iBeginLog,iEndLog,100);
                                end

                                logVideoCropped= cropLog(logVideo,iBeginVideo,iEndVideo,100);

                                logMerged     = appendStructs(logMerged,logVideoCropped);

                                % Create .txt file for YOLO groundTruth
                                if exportTXTFile && ~isempty(contains(videoFCNames,timeCharArray{i}(1:end-2)))
                                    videoFileName = videoFCNames{contains(videoFCNames,timeCharArray{i}(1:end-2))};
                                    m2Y(logMerged,fullfile(currCanapePath,videoFileName));
                                end
                            else
                                fprintf('\n Empty Video File !\n');
                            end
                        else
                            fprintf('\n Empty Video File !\n');
                        end
                    end
                    % Load GPS Recorder
                    if ~isempty(gpsFiles) && any(contains(gpsFiles,timeCharArray(i,:)))
                        recorderGPSFile = gpsFiles{contains(gpsFiles,timeCharArray(i,:))};
                        fprintf('\t\t Converting GPS : %s',recorderGPSFile);
                        logGPS      = readMDF_DM_MD(fullfile(currCanapePath,recorderGPSFile), 1, 0.01, 0,{'GPS_x','GPS_y'});
                        fprintf('\t -> Done !\n');
            %             toc
                        minPosixTime = max(logGPS.t_posixTime(1),logMerged.t_posixTime(1));
                        maxPosixTime = min(logGPS.t_posixTime(end),logMerged.t_posixTime(end));

                        [closestSample iBeginGPS] = min(abs(logGPS.t_posixTime-minPosixTime));
                        [closestSample iEndGPS]   = min(abs(logGPS.t_posixTime-maxPosixTime));
                        [closestSample iBeginLog] = min(abs(logMerged.t_posixTime-minPosixTime));
                        [closestSample iEndLog]   = min(abs(logMerged.t_posixTime-maxPosixTime));
                        if iEndGPS-iBeginGPS ~= iEndLog-iBeginLog
                            if iEndGPS-iBeginGPS == iEndLog-iBeginLog+1
                                iEndGPS = iEndGPS-1;
                            elseif iEndGPS-iBeginGPS == iEndLog-iBeginLog-1
                                iBeginGPS = iBeginGPS-1;
                            else
                                error('GPS log and Log have more than sample time delay !');
                            end
                        end
                        if iBeginLog>1 || iEndLog<size(logMerged.t,1)
                            logMerged = cropLog(logMerged,iBeginLog,iEndLog,100);
                        end

                        logGPSCropped= cropLog(logGPS,iBeginGPS,iEndGPS,100);

                        logMerged     = appendStructs(logMerged,logGPSCropped);
                    end
                    % save the log
                    fprintf('\t Saving %s.mat...',timeCharArray{i,3});
                    save(fullfile(currCanapePath,[timeCharArray{i,3} '.mat']),'-struct','logMerged');
                    fprintf('\t -> Done !\n');
        %             toc
                end
            end
            if concatenateLogsAnswer == 1
                concatenateLogs(currCanapePath,20);
            end
            end
        %     msgbox('Conversion finished !','Done');

        %% MergeVideos
        files = dir(fullfile(subFolderPath, '*.avi'));
        fileNames = {files.name};
        
        % Input file name
        fcFile = contains(fileNames, '_FC.avi');
        rcFile = contains(fileNames, '_RC.avi');
        
        if any(fcFile) && any(rcFile)
            file1 = fileNames{find(fcFile, 1)};
            file2 = fileNames{find(rcFile, 1)};
            
            video1 = VideoReader(fullfile(subFolderPath, file1));
            video2 = VideoReader(fullfile(subFolderPath, file2));
            
            % Output file name
            outputBaseName = regexprep(file1, 'FC\.avi$', 'CONTEXT.avi');
            outputFileName = fullfile(subFolderPath, outputBaseName);
            
            if exist(outputFileName, 'file')
                warning(['File ', outputFileName, ' already exists.']);
            else
                outputVideo = VideoWriter(outputFileName, 'Motion JPEG AVI');
                outputVideo.FrameRate = video1.FrameRate;
                open(outputVideo);
                
                % Set PIP size and position
                pipWidth = round(video2.Width / pi);
                pipHeight = round(video2.Height / pi);
                pipX = 10; 
                pipY = 10;
                
                totalPipFrames = video2.NumFrames;
                increment = video1.NumFrames / totalPipFrames;
                pipFrameIndex = 1;
                lastPipFrame = readFrame(video2);
                
                tic;
                lastPercentUpdate = 0;
                for k = 1:video1.NumFrames
                    frame1 = readFrame(video1);
                    if k >= round(pipFrameIndex)
                        if hasFrame(video2)
                            lastPipFrame = readFrame(video2);
                            pipFrameIndex = pipFrameIndex + increment;
                        end
                    end
                    
                    resizedPipFrame = imresize(lastPipFrame, [pipHeight pipWidth]);
                    frame1(pipY:(pipY + pipHeight - 1), pipX:(pipX + pipWidth - 1), :) = resizedPipFrame;
                    
                    writeVideo(outputVideo, frame1);
                    
                    currentProgress = round(100 * k / video1.NumFrames);
                    elapsedTime = toc;
                    framesPerSecond = k / elapsedTime;
                    estimatedRemainingTime = (video1.NumFrames - k) / framesPerSecond;
                    
%                     if currentProgress > lastPercentUpdate
                        waitbar(k / video1.NumFrames, h, sprintf('Rendering %s: %d%% - Est. Time: %.2f sec', subFolderNames{idx}, currentProgress, estimatedRemainingTime));
                        lastPercentUpdate = currentProgress;
%                     end
                end
                
                close(outputVideo);
                disp(['VIDEO MERGING COMPLETE: ', subFolderNames{idx}]);
            end
        else
            missingFiles = '';
            if ~any(fcFile)
                missingFiles = [missingFiles '_FC.avi '];
            end
            if ~any(rcFile)
                missingFiles = [missingFiles '_RC.avi'];
            end
            disp(['Missing files in ', subFolderNames{idx}, ': ', missingFiles]);
        end
    end
    
    close(h);
    
    %% local functions
    function Fichiers = filesearch(chemin, extension, avecSousDossiers)
if nargin < 2
    error('Pas assez d''arguments');
elseif nargin < 3
    avecSousDossiers = 1; % On séléctionne aussi tous les sous-dossiers
else
    if avecSousDossiers~=0 && avecSousDossiers~=1
        error('La variable ''avecSousDossiers'' doit être égale à 0 ou 1 !');
    end
end
FichiersDossierNiveau0 = dir(fullfile(chemin,'*'));
if ~numel(FichiersDossierNiveau0)
    error(['Le dossier' chemin 'est vide.']);
end

if avecSousDossiers == 0
    dirFlags    = [FichiersDossierNiveau0.isdir] &...
                   ~strcmp({FichiersDossierNiveau0.name},'.') & ~strcmp({FichiersDossierNiveau0.name},'..');
    FichiersDossierNiveau0 = FichiersDossierNiveau0(~dirFlags);
end

NbFichiers = 0;
Fichiers = struct('path', {}, 'name', {}, 'date', {}, 'bytes', {});
for f0 = 3:numel(FichiersDossierNiveau0)
    switch FichiersDossierNiveau0(f0).isdir
        case 1
            FichiersDossierNiveau1 = dir(fullfile(chemin, FichiersDossierNiveau0(f0).name));
            for f1 = 3:numel(FichiersDossierNiveau1)
                switch FichiersDossierNiveau1(f1).isdir
                    case 1
                        FichiersDossierNiveau2 = dir(fullfile(chemin, FichiersDossierNiveau0(f0).name, FichiersDossierNiveau1(f1).name));
                        for f2 = 3:numel(FichiersDossierNiveau2)
                            switch FichiersDossierNiveau2(f2).isdir
                                case 1
                                case 0
                                    [~, ~, ext] = fileparts(FichiersDossierNiveau2(f2).name);
                                    if strcmpi(ext,['.' extension])
                                        NbFichiers = NbFichiers + 1;
                                        Fichiers(NbFichiers).name = FichiersDossierNiveau2(f2).name ;
                                        Fichiers(NbFichiers).path = fullfile(chemin , FichiersDossierNiveau0(f0).name , FichiersDossierNiveau1(f1).name);
                                        Fichiers(NbFichiers).date = FichiersDossierNiveau2(f2).date;
                                        Fichiers(NbFichiers).bytes = FichiersDossierNiveau2(f2).bytes;
                                    end
                            end
                        end
                    case 0
                        [~, ~, ext] = fileparts(FichiersDossierNiveau1(f1).name);
                        if strcmpi(ext,['.' extension])
                            NbFichiers = NbFichiers + 1;
                            Fichiers(NbFichiers).name = FichiersDossierNiveau1(f1).name ;
                            Fichiers(NbFichiers).path = fullfile( chemin , FichiersDossierNiveau0(f0).name );
                            Fichiers(NbFichiers).date = FichiersDossierNiveau1(f1).date;
                            Fichiers(NbFichiers).bytes = FichiersDossierNiveau1(f1).bytes;
                        end
                    end
                end
            case 0
                [~, ~, ext] = fileparts(FichiersDossierNiveau0(f0).name);
                if strcmpi(ext,['.' extension])
                    NbFichiers = NbFichiers + 1;
                    Fichiers(NbFichiers).name = FichiersDossierNiveau0(f0).name;
                    Fichiers(NbFichiers).path = chemin;
                    Fichiers(NbFichiers).date = FichiersDossierNiveau0(f0).date;
                    Fichiers(NbFichiers).bytes = FichiersDossierNiveau0(f0).bytes;
                end
        end
    end
    if ~numel(NbFichiers)
        error(['Les dossiers ne contiennent aucun fichier *.' extension]);
    end

end

function log = readMDF_DM_MD(fileFullName,Interpolation, SampleTime, SetStartTime0,listVars,listNames)
    
    if nargin < 6 % listNames not defined
        listNames = {};
    end
    if nargin < 5 % listVars not defined
        listVars = {};
    end
    
    [finalDatas,  FieldMatrix] = Read_MDF_DM(fileFullName, Interpolation, SampleTime, SetStartTime0,listVars,listNames);
    if ~isempty(finalDatas)
        log = finishMDFConversion(finalDatas);
        log = compuCoefs(log,'RadarLongitudinalDistanceObject',0.078125,1);
        log = compuCoefs(log,'RadarLateralDistanceObject',0.078125,1);
        log = compuCoefs(log,'RadarAbsoluteLongiVelocityObject',0.00390625,1);
        log = compuCoefs(log,'RadarAbsoluteLatVelocityObject',0.00390625,1);
        log = compuCoefs(log,'RadarRelativeLatVelocityObject',0.00390625,1);
        log = compuCoefs(log,'RadarAbsoluteLongiAccelerationObject',0.00048828125,1);
        log = compuCoefs(log,'RadarAbsoluteLatAccelerationObject',0.00048828125,1);
        log = compuCoefs(log,'RadarRelativeLatAccelerationObject',0.00048828125,1);
        log = compuCoefs(log,'RadarLengthObject',0.0078125,1);
        log = compuCoefs(log,'RadarWidthObject',0.0078125,1);
        log = compuCoefs(log,'RadarHeadingAngleObject',0.00012207,1);
    else
        log = struct();
    end
end


function [finalDatas,  FieldMatrix] = Read_MDF_DM(curPath, Interpolation, SampleTime, SetStartTime0,listVars,listNames)
if nargin < 6 % listNames not defined
    listNames = {};
end
if nargin < 5 % listVars not defined
    listVars = {};
end

listStatsWith = cellfun(@(x) x(1:end-1),listVars(endsWith(listVars,'*')),'UniformOutput',false);
listEndsWith  = cellfun(@(x) x(2:end),listVars(startsWith(listVars,'*')),'UniformOutput',false);
listStartAndEndsWith = cellfun(@(x) x(2:end-1),listVars(startsWith(listVars,'*')&endsWith(listVars,'*')),'UniformOutput',false);
listMatch     = listVars(~contains(listVars,'*'));


FieldMatrix = {};
finalDatas = {};
% Read_MDF

% [FileName,PathName,~] = uigetfile('*.mdf;*.mf4','MultiSelect','on');
warning('off','all');
curBarH = waitbar(0,'Parsing MDF...');
curBarHb=findobj(curBarH,'Type','figure');
curBarHt = get(get(curBarHb,'currentaxes'),'title');
try
    mdfObj = mdf(curPath);
catch ME
    warning(ME.message);
    fprintf('\t --> Ignored.\n');
    close(curBarH);
    return;
end

% Parse all available rasters
curNames = get(mdfObj, 'ChannelNames');
curSizes = cellfun(@length, curNames);
timeVectors = cell(length(curNames(:,1)), 2);
timeVectorsErrors = zeros(length(curNames(:,1)), 1);
dataTable = cell(sum(curSizes), 3);
dataTableErrors = zeros(sum(curSizes), 1);
last_index = 1;
waitbar(0.1,curBarH);

%% Modified for signal name more than 63 char  
for i=1:length(curNames(:,1))
    timeVectors{i, 1} = ['t' num2str(i)];
    try timeVectors{i, 2} = read(mdfObj, i, mdfObj.ChannelNames{i}{curSizes(i)}, 'OutputFormat','Vector');
        if isempty(timeVectors{i,2})
            timeVectorsErrors(i, 1) = 1;
        end
    catch
        timeVectorsErrors(i, 1) = 1;
    end
 
    for j=1:(curSizes(i) - 1)
        curSignalName = curNames{i}{j};
        if length(curSignalName) > 63 % Apply the method of 'Read_MDF_DM_AC.m' 
            if contains(curSignalName,'.')
                parts = strsplit(curSignalName,'.');
                parts = strrep(strrep(parts,']',''),'[','_');
                if (length(char(parts(end))) + length(char(parts(end-1))) < 62)
                    curSignalName = [char(parts(end-1)) '_' char(parts(end))];
                else
                    curSignalName = char(parts(end));
                end
            else
                curSignalName = curSignalName(end-62:end);
            end
        end
        
        curNames{i}{j} = strrep(curSignalName, '[', '_');  
        curNames{i}{j} = strrep(curNames{i}{j}, ']', '');  

        if isempty(listVars) ||  wildcard(strrep(curNames{i}{j},'.','_'),listStatsWith,listEndsWith,listStartAndEndsWith,listMatch)
            if (Interpolation)
                dataTable{last_index, 1} = curNames{i}{j};
            else
                dataTable{last_index, 1} = [curNames{i}{j} '_' timeVectors{i, 1}];
            end
            dataTable{last_index, 3} = timeVectors{i,2};
            try dataTable{last_index, 2} = read(mdfObj, i, mdfObj.ChannelNames{i}{j}, 'OutputFormat','Vector');
            catch
                dataTableErrors(last_index, 1) = 1;
            end
            last_index = last_index + 1;
        end
    end
    waitbar(0.1+i/length(curNames(:,1))*0.8,curBarH);
end

dataTable = dataTable(~dataTableErrors, :);
dataTable = dataTable(~cellfun(@isempty, dataTable(:, 1)), :);
timeVectors = timeVectors(~timeVectorsErrors, :);
timeVectors = timeVectors(~cellfun(@isempty, timeVectors(:, 1)), :);

curTimes2Print = cellfun(@(x,y) [x,'_Length_of_t=',num2str(length(y))],timeVectors(:,1),timeVectors(:,2),'UniformOutput',false);

if (~Interpolation)
    [ChosenTimeIdx,~] = listdlg('ListString',curTimes2Print,'Name','Temps','ListSize',[300 400]);
    t=timeVectors{ChosenTimeIdx,2}(:,1);
    dataTable = dataTable(:,1:2);
    [~,idx] = sort(upper(dataTable(:,1)));
    finalDatas = [dataTable(idx, :) ; timeVectors ; {'t', t}];
else
    curBarHt.String = 'Interpolation...';
    ValidIdx = cellfun(@(x, y) ~isempty(x) && ~isempty(y) && length(x) == length(y), dataTable(:,2) , dataTable(:,3));
    MaxT = max(cellfun(@(x) x(end), timeVectors(:,2)));
    MinT = min(cellfun(@(x) x(1), timeVectors(:,2)));
    NewTime = (MinT:SampleTime:MaxT)';
    dataTable = dataTable(ValidIdx, :);
    %     dataTable(:,3) = cellfun(@(x) x-x(1), dataTable(:,3), 'UniformOutput', false);
    idxUnitary = cellfun(@(x) length(x) == 1,dataTable(:,2));
    idxCell    = cellfun('isclass',dataTable(:,2),'cell');
    if (~isempty(find(idxUnitary,1)))
        unitTable = dataTable(idxUnitary & ~idxCell, :);
        unitTable(:,3) = cellfun(@(x) NewTime, unitTable(:,3), 'UniformOutput', false);
        unitTable(:,2) = cellfun(@(x) x*ones(length(NewTime), 1), unitTable(:,2), 'UniformOutput', false);
    end
    
    dataTable = dataTable(~idxUnitary & ~idxCell, :);
%     dataTable(:,2) = cellfun(@(x, y) interp1(double(x), double(y), NewTime, 'linear',y(1)),dataTable(:,3), dataTable(:,2), 'UniformOutput', false);
    dataTable(:,2) = cellfun(@(x,y) interpCellFun(x,y,NewTime),dataTable(:,3), dataTable(:,2), 'UniformOutput', false);
    if (~isempty(find(idxUnitary & ~idxCell,1)))
        dataTable = [dataTable ; unitTable];
    end
    if (SetStartTime0)
        t = NewTime - NewTime(1);
    else
        t = NewTime;
    end
    dataTable = dataTable(:,1:2);
    [~,idx] = sort(upper(dataTable(:,1)));
    finalDatas = [dataTable(idx, :) ; {'t', t}];
end

%% adding for posixtime contained as initialTimestamp of mdfobj
% finalDatas = [finalDatas;{???????'t_posixtime',t+posixtime(mdfObj.InitialTimestamp)}???????];
if Interpolation
    finalDatas = [finalDatas;{'t_posixTime',posixtime(mdfObj.InitialTimestamp)+t}];
end

waitbar(1,curBarH);
close(curBarH)

warning('on','all');

% finalDatas(:,3) = cell(length(finalDatas(:,1)), 1);
if ~isempty(listNames)
    for ii=1:size(finalDatas,1)
        iFirstMatch = 0;
        nameFound = false;
        while ~nameFound && iFirstMatch<size(listNames,1)
            iFirstMatch = iFirstMatch+1;
            nameFound = contains(finalDatas{ii,1},listNames{iFirstMatch});
        end
        if nameFound
            finalDatas{ii,1} = finalDatas{ii,1}(regexp(finalDatas{ii,1},listNames{iFirstMatch},'matchcase','once'):end);
        end
    end
end
finalDatas(:,1) = matlab.lang.makeValidName(finalDatas(:,1));
finalDatas(:,3) = num2cell(ones(length(finalDatas(:,1)) , 1));
end
% FUNCTIONS

% Wildcar function -> filter signals according to listvar
function outputBool = wildcard(sigList,listStatsWith,listEndsWith,listStartAndEndsWith,listMatch)
    startWithBool        = any(cellfun(@(x) any(startsWith(sigList,x)),listStatsWith));
    endWithBool          = any(cellfun(@(x) any(endsWith(sigList,x)),listEndsWith));
    startAndEndsWithBool = any(cellfun(@(x) any(startsWith(sigList,x))& any(endsWith(sigList,x)),listStartAndEndsWith));
    listMatch            = any(cellfun(@(x) any(isequal(sigList,x)),listMatch));
    
    outputBool = startWithBool || endWithBool || startAndEndsWithBool || listMatch;
end

function yInterp = interpCellFun(x,y,NewTime)
    [xUnique iUnique] = unique(x);
%     offsetTime = xUnique(1)-NewTime(1);

    X = double(x(iUnique));%-offsetTime); % Current time, shiffted to 0
    V = double(y(iUnique)); % Current signal values
    Xq = NewTime; % Common interpolation time
    
    
    yInterp = interp1(X,V,Xq,'nearest','extrap');
    interpRange = Xq > min(X) & Xq < max(X);
    yInterp(interpRange) = interp1(X,V,Xq(interpRange),'next');
end

function log = finishMDFConversion(dataMatrix)
    log = struct();
    for i=1:size(dataMatrix,1)
        log.(dataMatrix{i,1}) = dataMatrix{i,2};
    end
end

function log = compuCoefs(log,sigLabel,compuNum,compuDen)
    % Get all log fields
    fields = fieldnames(log);
    % Find all signals containing specified label
    correspSig = fields(find(contains(fields,sigLabel)));
    
    for i = 1:size(correspSig,1)
        log.(correspSig{i})     = log.(correspSig{i})/(compuNum/compuDen);
    end
end

function log = cropLog(log,iBegin,iEnd,fCan)
    
    % get all signal names
    fields = fieldnames(log);
    for f=1:length(fields)
        if isequal(fields{f},'t')
            log.(fields{f})    = [0:(iEnd-iBegin)]'./fCan;
        elseif size(log.(fields{f}),1)>=iEnd-iBegin+1
            log.(fields{f})    = log.(fields{f})(iBegin:iEnd);
        else
            log.(fields{f})    = double(log.(fields{f}));
        end
    end
end

function s3 = appendStructs(s1,s2,prefix)
    if nargin < 3
        prefix = '';
    end
    s3 = s1;
    
    f2 = fieldnames(s2);
    for i=1:size(f2,1)
        s3.([prefix f2{i}]) = s2.(f2{i});
    end
end
```

