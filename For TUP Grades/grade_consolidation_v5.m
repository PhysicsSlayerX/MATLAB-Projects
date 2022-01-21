close all;
clc;
clear all;

% START TIMER
time1 = clock; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters
charPos         = 12;
charAdj         = 2;
numberOfQuizzes = 7;
namePrefix      = "TQ";

% Raw Data Configuration
nameColIndex    = 8;
scoreColIndex   = 6;
sectionColIndex = 11;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMPORT ALL THE NECESSARY FILES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fileNames = strings(numberOfQuizzes,1);
for i = 1:numberOfQuizzes
    fileNames(i,:) = strcat(namePrefix, num2str(i)); 
end


% Create a Struct File for storing imported table data
f0 = 'ThermoQuiz';      v0 = zeros(1,size(fileNames,1));

MainTable = struct(f0,v0);
for m = 1:size(fileNames,1)
MainTable(m).ThermoQuiz   = readtable(fileNames(m), 'ReadVariableNames', false);
end

% Create a Struct File for storing exported table data
a0 = 'ExportedData';      b0 = zeros(1,size(fileNames,1));

ExportedTable = struct(a0, b0);

% Import the Alphabetical List of Names
NamesList   = readtable('NamesThermoAlphabetical.csv', 'ReadVariableNames', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS LOOP FOR ALL THE FILES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for N = 1:size(fileNames,1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FILE PROCESSING
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Process the files
    MainFile          = MainTable(N).ThermoQuiz;
    MF_Cell           = table2cell(MainFile);
    ColumnsOfInterest = MF_Cell(:,[nameColIndex scoreColIndex sectionColIndex]);
    NamesList_c       = table2cell(NamesList);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DATA CLEANING AND FORMATTING
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Check and Remove Duplicates
    withDuplicates      = ColumnsOfInterest;
    nameCol             = convertCharsToStrings(withDuplicates(:,1));
    scoreCol            = num2str(double(cell2mat(withDuplicates(:,2))));
    sectionCol          = convertCharsToStrings(withDuplicates(:,3));
    [~, index]          = unique(strcat(nameCol,scoreCol ,sectionCol) , 'rows');
    index               = sort(index);
    withoutDuplicates   = withDuplicates(index,1:3);

    % Remove the spaces and commas in the names
    withoutDuplicates(:,1) = replace(withoutDuplicates(:,1), ...
                            {'.' ',' ' '},...
                            {'','',''});
    NamesList_c(:,1) = replace(NamesList_c(:,1), ...
                            {'.' ',' ' '},...
                            {'','',''});


    % Sorting Data According to Section and Names (Soon)


    % Determine Loop Sizes
    S = size(NamesList,1);
    D = size(withoutDuplicates,1);

    % Display Preliminary Data Infomation
    diffTakers = S - D;
    disp(['The total number of students: ' num2str(S)]);
    disp(['The total number of takers: ' num2str(D)]);
    disp(['The difference: ' num2str(diffTakers)]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PRE-ALLOCATION 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Pre-allocate arrays
    truth_table = zeros(S, S);
    matched_scores = zeros(S, 1);
    names = strings(S, 1);
    section_label = strings(S, 1);

    % Matching the names and scores
    for i = 1:S

        a = cell2mat(NamesList_c(i,1));
        counter = 0;

        for j = 1:D

            b = cell2mat(withoutDuplicates(j,1));
            c = double(cell2mat(withoutDuplicates(j,2)));
            d = convertCharsToStrings(cell2mat(withoutDuplicates(j,3)));


            if any(isnan(c)) || any(isempty(b)) || any(isempty(a))
                break;
            end

            %Create a regular expression
            %This expression matches any character except a whitespace, comma, period, semicolon, or colon 
            exp = '[^\s.,;:]*';
            %Find the matches within the string
            if length(b) < charPos || length(a) < charPos
                b1 = regexpi(b(1:charPos - charAdj), exp, 'match');
                %Concatenate all matches into a single string
                b1 = [b1{:}];
                %Repeat above for the other string
                a1 = regexpi(a(1:charPos - charAdj), exp, 'match');
                a1 = [a1{:}];
            else
                b1 = regexpi(b(1:charPos), exp, 'match');
                %Concatenate all matches into a single string
                b1 = [b1{:}];
                %Repeat above for the other string
                a1 = regexpi(a(1:charPos), exp, 'match');
                a1 = [a1{:}];
            end
            
            %Compare the modified characters
            truth_table(i,j) = strcmpi(a1, b1);
            disp([a1 ,"  ", b1])
            if truth_table(i,j)
               matched_scores(i,:) = c;
               names(i,:)          = convertCharsToStrings(a);
               section_label(i,:)  = d;
            end
            counter = counter + double(truth_table(i,j));

        end
        disp(counter);
        if counter < 1
            matched_scores(i,:) = 0;
            names(i,:)          = convertCharsToStrings(a);
            section_label(i,:)  = d;
        end
    end
    truth_table = sparse(truth_table);
    ExportedTable(N).ExportedData = horzcat(names, matched_scores, section_label);
    %disp(truth_table);
end


dataMatrix      = cell(size(NamesList, 1), size(fileNames, 1)+2);
dataMatrix(:,1) = table2cell(NamesList(:,1));
dataMatrix(:,2) = cellstr(ExportedTable(1).ExportedData(:,3));

% Post-processing

% First two columns are list of names and the section. Start at index 3
for N = 3:size(fileNames, 1)+2
    dataMatrix(:, N) = cellstr(ExportedTable(N-2).ExportedData(:,2));
end

writecell(dataMatrix,'Scores_TQ_v2.csv') 


% STOP TIMER
time2 = clock;
t = etime(time2,time1);
disp(['Elapsed time is ' num2str(t) ' seconds.']);
disp(['Elapsed time is ' num2str(t/60) ' minutes.']);
disp(['Elapsed time is ' num2str(t/60/60) ' hours.']);
disp(['Elapsed time is ' num2str(t/60/60/24) ' days.']);
