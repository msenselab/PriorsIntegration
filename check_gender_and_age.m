% Check the participants' gender information. Return a list of genders in
% certain folder. Just click into data folder. 

%% get gender

files=dir('.\*.mat')
g=[];
for i=1:length(files)
    load(files(i).name);
%     g{i}=a.sGender;
    g{i}=expInfo.sex;
end
g


% check how many males 
nm=0;
for i=1:length(g)
    if g{i}=='M'
        nm=nm+1;
    end
end
nm

%% get age

files=dir('.\*.mat')
a=[];
for i=1:length(files)
    load(files(i).name);
    a(i)=expInfo.age;
end
a

% calculate mean age and standard deviation
mean(a)
std(a)

