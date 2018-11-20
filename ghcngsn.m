%%%% 
%% before starting, you'll need to download the data from: 
%% then, unzip them, and go to the directory gsn/ghcnd_gsn/
%% do this on linux or make to make a text file with the names of all the data files in it: 
%% ls -1 > fnames.txt

% also change this to your directory
directory='/path/to/files/'


%%% from here it might just run. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ghcninv=importdata([directory 'ghcnd-inventory.txt']);
%         data: [nstationx2 double]
%     textdata: {nstationx4 cell}

% FORMAT OF "ghcnd-inventory.txt"
% 
% ------------------------------
% Variable   Columns   Type
% ------------------------------
% ID            1-11   Character
% LATITUDE     13-20   Real
% LONGITUDE    22-30   Real
% ELEMENT      32-35   Character
% FIRSTYEAR    37-40   Integer
% LASTYEAR     42-45   Integer
% ------------------------------

% find stations with precip data 
goodones=find(strcmp('PRCP',ghcninv.textdata(:,4)));
years=ghcninv.data(goodones,:);

invids=ghcninv.textdata(goodones,1);
statlats=ghcninv.textdata(goodones,2);
statlons=ghcninv.textdata(goodones,3);

for i=1:length(goodones)
    lats(i)=str2num(statlats{i});
    lons(i)=str2num(statlons{i});
end

figure(1);clf
hold on 
plot(lons,lats,' kx')
plot(lons,lats,' rx')


dateres=length(invids) 


ghcnstations=importdata([directory 'ghcnd-stations.txt']);

% FORMAT OF "ghcnd-stations.txt" 
% ------------------------------
% Variable   Columns   Type
% ------------------------------
% ID            1-11   Character
% LATITUDE     13-20   Real
% LONGITUDE    22-30   Real
% ELEVATION    32-37   Real
% STATE        39-40   Character
% NAME         42-71   Character
% GSN FLAG     73-75   Character
% HCN/CRN FLAG 77-79   Character
% WMO ID       81-85   Character
% ------------------------------

for i=1:size(ghcnstations.textdata,1)
    t=textscan(ghcnstations.textdata{i},'%s %f %f %f %31c %s %s');
    sid{i}=t{1};
    elev{i}=t{4};
    elev{i}=t{4};
    name{i}=t{5};
    gsn{i}=t{6};
    theend{i}=t{7};
end

for i=1:size(ghcnstations.textdata,1)
    t=textscan(ghcnstations.textdata{i},'%s %f %f %f %31c %s %s');
end



count=0;
for i=1:length(theend)
    if ~isempty(theend{i}); 
        count=count+1;
        notempty(count)=i;
    end
end

count=0;
for i=1:length(theend)
    if strcmp('GSN',gsn{i});
        count=count+1;
        isgsn(count)=i;
    end
end
isgsn=isgsn(:);

for i=1:length(notempty)
    if strcmp('GSN',theend{notempty(i)});
        isgsn(end+1)=notempty(i);
    end
end

isgsn=sort(isgsn);
sid=sid(isgsn);
elev=elev(isgsn);

count=0;
for i=1:length(isgsn)    
    in=find(strcmp(sid{i},invids));
    if ~isempty(in)
        count=count+1;
        all2gsn(count)=in;
    end
end

gsnlats=lats(all2gsn);
gsnlons=lons(all2gsn);
gsnids=invids(all2gsn);

plot(lons,lats,' kx',gsnlons,gsnlats,' rx')
legend('all ghcnd','gsn','location','northoutside')
xlim([-180 180])
ylim([-90 90])
set(gcf,'units','centimeters','paperpositionmode','auto');
set(gcf,'position',[10 14 8.4 8.4]); %% one col: 8.4 cm. two col 16.9.  max height 23.7.


% check to make sure we really have all these files 
fnames=importdata([directory '/gsn/ghcnd_gsn/fnames.txt']);

for i=1:(size(fnames,1)-1)
    t=textscan(fnames{i},'%11c %s');
    fileids{i}=t{1};
end

count=0;
missingfile=[]
for i=1:length(gsnids)
    if isempty(find(strcmp(gsnids(i),fileids(:))))
        count=count+1;
        missingfile(count)=i;
    end
end
missingfile


gsnids(missingfile)=[];
gsnlats(missingfile)=[];
gsnlons(missingfile)=[];
all2gsn(missingfile)=[];

years=ghcninv.data(goodones,:);
gsnstartyears=years(all2gsn,1);
gsnendyears=years(all2gsn,2);
firststartyear=min(gsnstartyears)
laststartyear=max(gsnstartyears)

startyears=firststartyear:laststartyear;
for i=1:length(startyears)
    nstaty(i)=sum(gsnstartyears<=startyears(i));
end

figure(2);
plot(startyears,nstaty,'ok-')


startyear=1893; % this year there's a big jump in the number of
                % stations online 
endyear=str2num(datestr(now,10))-1;

prcp=NaN([length(gsnids) length(startyear:endyear+1) 12 31]);

for stati=1:length(gsnids)
    disp([num2str(stati) ' ' gsnids{stati}])
    fid=fopen(strcat(directory,'gsn/ghcnd_gsn/',gsnids{stati},'.dly'));
    tline=fgetl(fid);
    while ischar(tline)
        year=str2num(tline(12:15));
        if year>(startyear-1)
            month=str2num(tline(16:17));
            variable=tline(18:21);
            switch variable
                case 'PRCP' 
                    prcp(stati,year-(startyear-1),month,:)=readline_ghcnd_prcp(tline);
            end
        end
        tline=fgetl(fid);
    end
    fclose(fid);
end

prcpraw=prcp; 


pnz=prcp;
pnz(prcp>0)=1;
prcp(prcp==eps)=0;



years=startyear:endyear+1;

save gsndailydata.mat prcp pnz gsnids gsnlats gsnlons ...
    years gsnstartyears gsnendyears

                                                                      

                                                                                                                   


