close all;
clc;
clear all;

% START TIMER
time1 = clock; 

% Read the table/s
T = readtable('SampleCSV.csv');
V = readtable('NamesThermoAlphabetical.csv');
MainFile = readtable("TQ1");

% Create a Struct File for storing table data
f0 = 'ThermoQuiz';      v0 = zeros(500,500);
MainTable = struct(f0,v0);

% Sorting Data According to Section and Names

% Convert Table to Cell to Array
wd = table2cell(T);
col2 = convertCharsToStrings(wd(:,2));
col3 = num2str(double(cell2mat(wd(:,3))));
col4 = convertCharsToStrings(wd(:,4));
[~, index] = unique(strcat(col2,col3,col4) , 'rows');
index = sort(index);
noDup = wd(index,2:4);
S = size(V,1);
D = size(noDup,1);

% Pre-allocate arrays
truth_table = zeros(S, S);
matched_scores = zeros(S, 1);
names = strings(S, 1);
section_label = strings(S, 1);
% Character Splicing Adjustments
charPos = 12;
charAdj = 1;

% Matching the names and scores
for i = 1:S
    a = cell2mat(V{i,1});
    counter = 0;
    
    for j = 1:D
        
        b = cell2mat(noDup(j,1));
        c = double(cell2mat(noDup(j,2)));
        d = cell2mat(noDup(j,3));
        
        if any(isnan(c)) || any(isempty(b)) || any(isempty(a))
            break;
        end

        %a = 'test';
        %b = 'te s.t';

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
           section_label(i,:)  = convertCharsToStrings(d);
        end
        counter = counter + double(truth_table(i,j));
        
    end
    disp(counter);
    if counter > 2
        matched_scores(i,:) = 0;
        names(i,:)          = convertCharsToStrings(a);
        section_label(i,:)  = convertCharsToStrings(d);
    end
end

% Post-processing
truth_table = sparse(truth_table);
Score_tally = horzcat(names, matched_scores, section_label);
disp(truth_table);

writematrix(Score_tally,'Scores_ThermoQuiz7.csv') 


% STOP TIMER
time2 = clock;
t = etime(time2,time1);
disp(['Elapsed time is ' num2str(t) ' seconds.']);
disp(['Elapsed time is ' num2str(t/60) ' minutes.']);
disp(['Elapsed time is ' num2str(t/60/60) ' hours.']);
disp(['Elapsed time is ' num2str(t/60/60/24) ' days.']);
