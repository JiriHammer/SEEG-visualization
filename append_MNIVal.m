function [ MNI, Vals ] = append_MNIVal( MNIAll, ValsAll )
%APPEND_MNIVAL udela jednu strukturu a pole z cells 
%   to se pak da pouzit do Jirkovy funkce main_brainPlot
%   jednotlive prvky MNIAll a ValsAll se ziskaji z CiEEGData.IntervalyResp (s odebranim posledniho rozmeru)
%   (c) Kamil Vlcek 28.3.2017

assert(numel(MNIAll)==numel(ValsAll), 'MNIAll a ValsAll musi mit stejny pocet prvku');
%nejdriv spocitam celkovy pocet kanalu pro prealokaci
channelsMNI  = 0;
channelsVals = 0;
for j = 1:numel(MNIAll)
    channelsMNI = channelsMNI + numel(MNIAll{j});
    channelsVals = channelsVals + size(ValsAll{j},1);
end
assert(channelsMNI == channelsVals, 'celkovy pocet kanalu musi byt stejny v MNI a Vals');
MNI = repmat(struct('MNI_x',0,'MNI_y',0,'MNI_z',0),channelsMNI,1); 
Vals = zeros( channelsVals,size(ValsAll{1},2)); % channels vs intervaly casu
%pak vsechny budouci kanaly postupne plnim
CH = 1;
for j = 1:numel(MNIAll)
    mnitemp =  MNIAll{j}; %pro jednodussi zapis
    valstemp = ValsAll{j};
    for ch = 1:numel(mnitemp);       
       MNI(CH).MNI_x =mnitemp(ch).MNI_x; 
       MNI(CH).MNI_y =mnitemp(ch).MNI_y; 
       MNI(CH).MNI_z =mnitemp(ch).MNI_z;        
       Vals(CH,:) = valstemp(ch,:);
       CH = CH + 1;
    end
end
end

