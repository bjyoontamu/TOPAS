%% TOPAS algorithm| Network-based RNA structural alignment algorithm
%-------------------------------------------------------------------------
%% TOPAS is RNA structural alignment through topological networks.
%  Input:	
%       input_sequences	the input sequences in FASTA format
%      	input_prob1     the input probability of base pairing for sequence1
%      	input_prob2     the input probability of base pairing for sequence2
%      	input_prob3     the input probability of base alignment with format
%      	input_parm1     the input parameter alpha for structure similarity
%      	input_parm2     the input parameter beta for connected similarity
%                       (0<= alpha+ beta<= 1)
%
% Output:
%        output_seq1		the output aligned sequence 1
%        output_seq2		the output aligned sequence 2
%% Author: Aky
%% version 1.3 08.28.2018
%--------------------------------------------------------------------------


function [Aseq1,Aseq2]= TOPAS(Fseqs,Fbp1,Fbp2,Fap12, Alpha, Beta)


%% Code parameter
global FAData;
global Seq1 Seq2;
global MPb1 Len1 Len2;
global PS1 SPS1 PS2 SPS2;

assert(Alpha>=0, 'Parameter incorrect!');
assert(Beta>=0, 'Parameter incorrect!');
assert(Alpha+Beta<=1, 'Parameter incorrect!');

fprintf('Read FASTA %s\n',Fseqs);
seqread(Fseqs);

NIt= 50;      % Num Max Iteration
Toler= 1e-2;  % Iternation tolerance
Thred= 0.01;  % Threshold Prob struct


%% Read Dataset
display('Load Model...');

ME= zeros(Len1,Len2);
PS1= zeros(Len1);
PS2= zeros(Len2);

%%  Read Prob model
assert(exist(Fap12,'file')==2, 'File %s not found!', Fap12);
FPData= dlmread(Fap12);
assert(length(FPData(1,:))==3, 'Format %s incorrect!', Fap12);
ME(sub2ind(size(ME),FPData(:,1),FPData(:,2)))=FPData(:,3);
ME= ME/sum(ME(:));

assert(exist(Fbp1,'file')== 2, 'File %s not found!',Fbp1);
FPData= dlmread(Fbp1);
assert(length(FPData(1,:))==3, 'Format %s incorrect!', Fbp1);
PS1(sub2ind(size(PS1),FPData(:,1),FPData(:,2)))=FPData(:,3);
%% Threshold Filter1
PS1(find(PS1<Thred))= 0;
PS1(find(eye(Len1)))= 0;
PS1= PS1+PS1';
SPS1= sum(PS1,2);

assert(exist(Fbp2,'file')== 2, 'File %s not found!',Fbp2);
FPData= dlmread(Fbp2);
assert(length(FPData(1,:))==3, 'Format %s incorrect!', Fbp2);
PS2(sub2ind(size(PS2),FPData(:,1),FPData(:,2)))=FPData(:,3);
%% Threshold Filter2
PS2(find(PS2<Thred))= 0;
PS2(find(eye(Len2)))= 0;
PS2= PS2+PS2';
SPS2= sum(PS2,2);

%% Run PINET
display('TOPAS Align...');
MPb2= 1+rand(Len1,Len2);
MPb2= MPb2/sum(MPb2(:));

for k=1:NIt
  MPb1= MPb2;

  for (a=1:Len1)
    for (b=1:Len2)
      AR1= up_str1(a,b);
      AR2= up_str2(a,b);
      ARs= Beta* AR1+ Alpha*AR2;
      MPb2(a,b)= ARs+ (1-Alpha-Beta)*ME(a,b);
    end
  end
  
  MPb2= MPb2/sum(MPb2(:));
  err= abs(MPb2-MPb1);
  if (sum(err(:))<Toler)
    break;
  end
end

%% Network alignment optimization
[Aseq1, Aseq2] = Needleman(FAData(1).Sequence,FAData(2).Sequence, MPb2 );


%_end Function


%% Read sequences
function [output]= seqread(Fname1)
global FAData;
global Seq1 Seq2 Len1 Len2;

%  Read Fasta file
assert(exist(Fname1,'file')==2, 'File %s not found!', Fname1);
FAData= fastaread(Fname1);

Seq1= nt2int(FAData(1).Sequence);
Seq2= nt2int(FAData(2).Sequence);
Len1= length(Seq1);
Len2= length(Seq2);


%_end function


%% Structural alignment

function [AR1]= up_str1(Id1,Id2)
global MPb1 Len1 Len2;

AR1= 0;

if (Id1>1&& Id2>1)
  AR1= AR1+ MPb1(Id1-1,Id2-1);
end  

if (Id1<Len1&& Id2<Len2)
  AR1= AR1+ MPb1(Id1+1,Id2+1);
end

AR1= AR1/2;

%_end function



function [AR2]= up_str2(Id1,Id2)
global MPb1 Len1 Len2;
global SPS1 PS1 SPS2 PS2;

AR2= 0;

for (a=1:Len1)
  if (PS1(Id1,a))
    for (b=1:Len2)
      if (PS2(Id2,b))
        S2= PS1(Id1,a)*PS2(Id2,b)/SPS1(a)/SPS2(b);
        AR2= AR2+ S2*MPb1(a,b);
      end
    end
  end
end

%_end function



%% Needleman Optimization
function [aseq1,aseq2] = Needleman(seq1,seq2,PS)

%% Parameter
n= length(seq1);
m= length(seq2);
w= 0; % Weight Gap

gamma = zeros(n+1, m+1); 
back = cell(n+1, m+1);

%% Fill the gamma matrix
for i=2:n+1
  for j=2:m+1
    sij= PS(i-1,j-1);
    a = [gamma(i-1,j-1)+sij, gamma(i,j-1)+w, gamma(i-1,j)+w];
    [smax, imax]= max(a);
    gamma(i,j)= smax;
    switch imax
      case 1
        back{i,j} = [i-1, j-1];
      case 2
        back{i,j} = [i, j-1];
      case 3
        back{i,j} = [i-1, j];
    end
  end
end


%% Back tracing
i= n+1;
j= m+1;
buf= cell(1,max([n, m])-1);
buf{1}= [i, j];
k= 2;
while (i>1 && j>1)
  ij= back{i,j};
  buf{k} = ij;
  k= k + 1;
  i= ij(1);
  j= ij(2);
end

%% Back init
aseq1=[];
aseq2=[];
ij= buf{length(buf)};
iold= ij(1);
jold= ij(2);


for (i=1:iold-1)
  aseq1= strcat(aseq1,seq1(i));
  aseq2= strcat(aseq2,'-');
end

for(j=1:jold-1)
  aseq1= strcat(aseq1,'-');
  aseq2= strcat(aseq2,seq2(j));
end

%% Back alignment
for k=length(buf)-1:-1:1
  ij= buf{k};
  i= ij(1);
  j= ij(2);
  if (i==iold)
    c1= '-';
  else
    c1 = seq1(i-1);
  end
  if (j==jold)
    c2= '-';
  else
    c2 = seq2(j-1);
  end
  
  iold= i;
  jold= j;
  aseq1= strcat(aseq1,c1);
  aseq2= strcat(aseq2,c2);

end


%% Alignment Optimization
for (i=iold:n)
  aseq1= strcat(aseq1,seq1(i));
  aseq2= strcat(aseq2,'-');
end

for(j=jold:m)
  aseq1= strcat(aseq1,'-');
  aseq2= strcat(aseq2,seq2(j));
end



%_end function

