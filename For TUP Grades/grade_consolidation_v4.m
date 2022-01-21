close all;
clc;
clear all;

% START TIMER
time1 = clock; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Character Splicing Adjustments
charPos = 12;
charAdj = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILE PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Process the files
MainFile          = readtable("TQ7_2", 'ReadVariableNames', false);
MF_Cell           = table2cell(MainFile);
ColumnsOfInterest = MF_Cell(:,[8 6 11]);

% Create a Struct File for storing table data
f0 = 'ThermoQuiz';      v0 = zeros(500,500);
MainTable = struct(f0,v0);

% Import the Alphabetical List of Names
NamesList   = readtable('NamesThermoAlphabetical.csv', 'ReadVariableNames', false);
NamesList_c = table2cell(NamesList);

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


% Sorting Data According to Section and Names


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
    
    a = cell2mat(NamesList_c(i,:));
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
        %Compare the modified strings
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

% Post-processing
truth_table = sparse(truth_table);
Score_tally = horzcat(names, matched_scores, section_label);
disp(truth_table);

writematrix(Score_tally,'Scores_TQ7.csv') 


% STOP TIMER
time2 = clock;
t = etime(time2,time1);
disp(['Elapsed time is ' num2str(t) ' seconds.']);
disp(['Elapsed time is ' num2str(t/60) ' minutes.']);
disp(['Elapsed time is ' num2str(t/60/60) ' hours.']);
disp(['Elapsed time is ' num2str(t/60/60/24) ' days.']);
