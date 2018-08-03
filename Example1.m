% TOPAS algorithm| Network-based RNA structural alignment algorithm
% Author: Aky Chen
% version 1.3 2017.09.01
% All right reserved

function Example1


%- Code parameter
Dhead= './dataset';
List= dir(Dhead);
List= List(arrayfun(@(a) a.name(1),List)~='.');


%-- PINET paramter
Alpha= 0.5;
Beta= 0.48;

for (a=1:length(List))
  Fpath= strcat(Dhead,'/',List(a).name);
  Fseqs= strcat(Fpath,'/seqs2.fa');
  Fbp1= strcat(Fpath,'/bp.0');
  Fbp2= strcat(Fpath,'/bp.1');
  Fap12= strcat(Fpath,'/ap.0.1');
  
  [aseq1,aseq2]= TOPAS(Fseqs,Fbp1,Fbp2,Fap12, Alpha, Beta);
  
  fprintf('Aligned Seqs:\n');
  fprintf('%s\n',aseq1);
  fprintf('%s\n\n',aseq2);

 
end






%_end simulation

