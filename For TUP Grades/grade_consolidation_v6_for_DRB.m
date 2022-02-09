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
charAdj         = 7;

% Item Details
numberOfItems   = [8; 5; 2];
prefix          = ["DRBQ"; "DRBPS"; "DRBE"];

% Raw Data Configuration
ColIndex            = [8 6 11];
ColIndex_exams      = [9 6 12];
ColIndex_probsets   = [9 6 12];

% List of Alphabetical Names
alphabeticalNames = 'NamesDRBAlphabetical.csv';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMPORT ALL THE NECESSARY FILES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath('D:\Research\TUP Important Files\2021-2022 First Sem\Quizzes\Quizzes_PS_Exams_Repo_v2')  

itemName = cell(size(numberOfItems,1),1);

% Pre-allocate the data containers
for j = 1: size(numberOfItems,1)
    % itemName(1, 1) -> Quiz Container
    % itemName(2, 1) -> ProbSet Container
    % itemName(3, 1) -> Exams Container
    itemName{j, 1} = strings(numberOfItems(j),1);
end

for k = 1: size(numberOfItems,1)
    for i = 1: numberOfItems(k)
        itemName{k,1}(i,:) = strcat(prefix(k,1), num2str(i));
    end
end
 quizzes  = itemName{1,1};
 probSets = itemName{2,1};
 exams    = itemName{3,1};
 
% Create a Struct File for storing imported table data
f0 = 'DRBQuiz';      v0 = zeros(1,size(quizzes,1));
f1 = 'DRBProbSets';  v1 = zeros(1,size(probSets,1));
f2 = 'DRBExams';     v2 = zeros(1,size(exams,1));
MainTable = struct(f0,v0,f1,v1,f2,v2);

for m = 1:size(quizzes,1)
    MainTable(m).DRBQuiz         = readtable(quizzes(m), 'ReadVariableNames', false);
end

for m = 1:size(probSets,1)
    MainTable(m).DRBProbSets     = readtable(probSets(m), 'ReadVariableNames', false);
end

for m = 1:size(exams,1)
    MainTable(m).DRBExams        = readtable(exams(m), 'ReadVariableNames', false);
end

% Import the Alphabetical List of Names
NamesList   = readtable(alphabeticalNames, 'ReadVariableNames', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS LOOP FOR ALL THE FILES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for reqIndex = 1 : size(itemName,1)
    
    % Create a Struct File for storing exported table data
    a0 = 'ExportedData';      b0 = zeros(1,size(ColIndex,2));
    
    ExportedTable = struct(a0, b0);
    
    for N = 1:numberOfItems(reqIndex)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % FILE PROCESSING
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Process the files
        if reqIndex == 1
            MainFile          = MainTable(N).DRBQuiz;
            MF_Cell           = table2cell(MainFile);
            ColumnsOfInterest = MF_Cell(:,ColIndex);
        elseif reqIndex == 2
            MainFile          = MainTable(N).DRBProbSets;
            MF_Cell           = table2cell(MainFile);
            ColumnsOfInterest = MF_Cell(:,ColIndex_probsets);
        else
            MainFile          = MainTable(N).DRBExams;
            MF_Cell           = table2cell(MainFile);
            ColumnsOfInterest = MF_Cell(:,ColIndex_exams);
        end  
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
        
        %pre-allocate withoutDuplicate array
        withoutDuplicates   = cell(size(index,1), size(ColIndex,2));
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
            
            D = size(withoutDuplicates,1);
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
                counter = counter + double(truth_table(i,j));
                if truth_table(i,j)
                   matched_scores(i,:) = c;
                   names(i,:)          = convertCharsToStrings(a);
                   section_label(i,:)  = d;
                   withoutDuplicates(j,:) = [];
                   break;
                end
               

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
        disp(truth_table);
    end

    disp(["Done processing the matching for " prefix(reqIndex)]);
    
    dataMatrix      = cell(size(NamesList, 1), size(itemName{reqIndex,1}, 1)+2);
    dataMatrix(:,1) = table2cell(NamesList(:,1));
    dataMatrix(:,2) = cellstr(ExportedTable(1).ExportedData(:,3));

    % Post-processing
    % First two columns are list of names and the section. Start at index 3
    for N = 3:size(itemName{reqIndex,1}, 1)+2
        dataMatrix(:, N) = cellstr(ExportedTable(N-2).ExportedData(:,2));
    end
    
    formattedFileName = strcat('Scores_', prefix(reqIndex) ,'.csv');
    writecell(dataMatrix, formattedFileName);
end




% STOP TIMER
time2 = clock;
t = etime(time2,time1);
disp(['Elapsed time is ' num2str(t) ' seconds.']);
disp(['Elapsed time is ' num2str(t/60) ' minutes.']);
disp(['Elapsed time is ' num2str(t/60/60) ' hours.']);
disp(['Elapsed time is ' num2str(t/60/60/24) ' days.']);
