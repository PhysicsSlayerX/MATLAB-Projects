close all;
clc;
clear all;

% Read the table
T = readtable('SampleCSV.csv');
MainFile = readtable("ThermoMidterms.csv");

% Convert Table to Cell to Array
wd = table2cell(T);
col2 = convertCharsToStrings(wd(:,2));
col3 = num2str(double(cell2mat(wd(:,3))));
[~, index] = unique(strcat(col2,col3) , 'rows');
index = sort(index);
noDup = wd(index,2:3);
S = 42;
D = size(noDup,1);

% Pre-allocate arrays
truth_table = zeros(S, S);
matched_scores = zeros(S, 1);
names = strings(S, 1);

% Matching the names and scores
for i = 1:S
    a = cell2mat(T{i,1});
    counter = 0;
    
    for j = 1:D
        
        b = cell2mat(noDup(j,1));
        c = double(cell2mat(noDup(j,2)));
        
        if any(isnan(c))
            break;
        end

        %a = 'test';
        %b = 'te s.t';

        %Create a regular expression
        %This expression matches any character except a whitespace, comma, period, semicolon, or colon 
        exp = '[^\s.,;:]*';
        %Find the matches within the string
        b1 = regexpi(b(1:13), exp, 'match');
        %Concatenate all matches into a single string
        b1 = [b1{:}];
        %Repeat above for the other string
        a1 = regexpi(a(1:13), exp, 'match');
        a1 = [a1{:}];
        %Compare the modified strings
        truth_table(i,j) = strcmpi(a1, b1);
        disp([a1 ,"  ", b1])
        if truth_table(i,j)
           matched_scores(i,:) = c;
           names(i,:)          = convertCharsToStrings(a);
        end
        counter = counter + double(truth_table(i,j));
        
    end
    disp(counter);
    if counter > 2
        matched_scores(i,:) = 0;
        names(i,:)          = convertCharsToStrings(a);
    end
end

% Post-processing
truth_table = sparse(truth_table);
Score_tally = horzcat(names, matched_scores);
disp(truth_table);