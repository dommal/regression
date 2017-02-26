%% This script calculates maximum adjusted R^2 within a number of
% possibilities checking for a specified p-value. Data has to imported from
% xls. file 

%% OUTPUT: R_max - maximum adjusted R2 within all combinations 
%terms- name of variables in the combination corresponding to R_max
%p-val-p-values for the corresponging R2=R_max, first row is intercept,
%next raws correspond to variabes (terms)
%%
clear 
clc
%INPUT
Dir='C:\Users\Name\XXXXXXXXXXXXXXXXXXXXXXX\'; % specify Directory for your Excel file
numVar=3; % specify number of your variables to be combined later
pValue=0.055; % specify p value. The script will eliminate all combinations which include any p-value bigger than specified
Y = xlsread('FileName.xlsx','SheetName','A2:A100'); % RESPONSE VARIABLE, it has to be a vector! 
X=xlsread('FileName.xlsx','SheetName','B2:D100'); %PREDICTOR VARIABLES,  matrix,  number of columns=numVar , length X=length Y
[num, string, raw] = xlsread('FileName.xlsx','SheetName','B1:D1');% reads names of your variables (top of your table in excel)
%END OF INPUT
v=[1:numVar];
for m=1:numVar
   C=nchoosek(v, m); % creates a matrix with all possible comination of your variables
   A = nchoosek(numVar,m); % numer of possible combinations

for i=1:A
    col=C(i,:) ;
    x=[X(:,col)];
    stats = regstats(Y,x,'linear',{'adjrsquare' 'tstat'});
    R(i) = getfield(stats, 'adjrsquare');
    Getval(i)= getfield(stats, 'tstat');
end

RR{m}=R; % adjusted R2
clear R;
Pval{m} = Getval; %p-values including intercept
clear Getval;

AA=Pval{1,m};
Values{m}= [AA.pval]; %p-values including intercept in the first row
RandP{m}=[Values{1,m} ; RR{1,m}]; % combined R and p-val, R at the end
   
    for bb=1:A
        ABC=RandP{1,m}(2:end-1,bb);
        if all(ABC<=pValue)==1; % select only cases where p val <=0.055 and put them in a new matrix
           GG{bb}=[RandP{1,m}(end,bb)];
        end
    end
  
OUTPUT{m}=GG ;
GG(:)=[];
  
DD=OUTPUT{1,m};
empties = cellfun('isempty',DD);
DD(empties) = {NaN}; % replace empty cells with NaN


R_mat{m} = cell2mat(DD); % converts from cell to number
    
end
R_mat=R_mat';

for yy=1:numVar
    AAA(yy)=length(R_mat{yy,1});
end
ddd=max(AAA); 

for iuy=1:numVar
    leng=length(R_mat{iuy,1});
    for ee=leng+1:ddd;
   R_mat{iuy,1}(1,leng+1:ddd)=NaN; %all vectors have to have the same length, those which are shorter are extended with NaNs
   
    end
end

all_R2 = cell2mat(R_mat);% all R2 
R_max=max(all_R2);
R_max=max(R_max) % maximum adjusted R2 within all possibilities with pval <=0.055

[i,j] = find(all_R2 == R_max)% i-numer of terms used , j-combination number

all_comb=nchoosek(v, i);
terms_num=all_comb(j,:);
terms=string(terms_num) % variables crrespondng to the best R2 combination
p_val=Values{1,i}(:,j) %p-values for the corresponging R2=R_max, first row is intercept
%% only run if needed: This section checks adjusted R2 for a specified number of variables combined at a time. e.g x=6 will check max adj R2 for all combinations using 6 variables.
%INPUT
x=3 % number from 1:numVar
%END OF INPUT
R_max2=max(all_R2(x,:)) % R2 max for number of terms =x
[k,l] = find(all_R2 == R_max2)
all_comb2=nchoosek(v, k);
terms_num2=all_comb2(l,:);
terms2=string(terms_num2) % variables crrespondng to the best R2 combination
p_val2=Values{1,k}(:,l) %p-values for the corresponging R2=R_max, first row is intercept