clear
clc
tic
%% setting  parameters
n=16; % Coarse-scale patch size
s=3;  % Species number
num=1024; % Fine-scale patch number
T = 1500; % generations
H = [0 1 0;0 0 1;1 0 0];% competition matrix
cf=100; % 100 independent simulations
CV=zeros(cf,1); 
load('ran.mat')  % random network matrix
network=M;
load('rana.mat')  % Fine-scale patch  matrix
for i1=1:cf
    im=i1;
    lamida = 2^(-10) ; % dispersal probability D
    local_rho=zeros(T,s*n);
    load('initialspecies.mat');
    for rep = 1:T
        rep
        for t=1:num
            focous_p = zeros(n,1);    % selecting empty patch
            for ig=1:n   % occupation patch
                focous_p(ig,:) = randperm(num,1)+(ig-1)*num;  % focal patch
            end
            I_species(focous_p)=0;
            species=I_species;
            pd1 = rand(n,1);
            pd = ones(n,1);pd(pd1<lamida)=0;
            %%
            for sc =1:n
                if pd(sc)==1
                    network2=a{1,1}; % Dispersal network of Species 1
                    network3=a{5,1}; % Dispersal network of Species 2
                    network4=a{10,1}; % Dispersal network of Species 3
                    death=focous_p(sc)-(sc-1)*num;
                    %% occupancy empty patch of Species 1
                    competitotrs1=find(network2(:,death)==1);
                    S_competitotrs1=I_species(competitotrs1,sc);
                    S_competitotrs1(S_competitotrs1~=1)=[];
                    %% occupancy empty patch of Species 2
                    competitotrs2=find(network3(:,death)==1);
                    S_competitotrs2=I_species(competitotrs2,sc);
                    S_competitotrs2(S_competitotrs2~=2)=[];
                    %%  occupancy empty patch of Species 3
                    competitotrs3=find(network4(:,death)==1);
                    S_competitotrs3=I_species(competitotrs3,sc);
                    S_competitotrs3(S_competitotrs3~=3)=[];
                    S_competitotrs=[ S_competitotrs1; S_competitotrs2; S_competitotrs3];
                    LM=length(S_competitotrs);
                    if  LM==0
                        species;
                    elseif  LM==1
                        species(focous_p(sc))=S_competitotrs(1,1);
                    else
                        competitotr= S_competitotrs(randperm(LM,2)');
                        competitotr1=sort(competitotr);
                        HH=H(competitotr1(2,1),competitotr1(1,1));
                        panduan=rand(1)<HH;
                        species(focous_p(sc))=competitotr1(panduan+1,1);
                    end
                    %% selecting species within Fine-scale patch occupy empty patch
                else
                    %%
                    linju = find(network(:,sc)==1);
                    link_local_cell=I_species(:,linju);
                    link_local_cell( link_local_cell==0)=[];
                    L1=length( link_local_cell);
                    x=randperm(L1,2);
                    competitotrs_exs=link_local_cell(x);
                    HH2=H( competitotrs_exs(1), competitotrs_exs(2));
                    species(focous_p(sc))=  competitotrs_exs(2-HH2);
                end
                %% selecting species within Coarse-scale patch occupy empty patch
            end
            I_species=species;
        end
        for isn = 1:n
            for isp=1:s
                local_rho(rep,isp+s*(isn-1)) =length(find(I_species(:,isn)==isp))/num;
            end
        end
    end
    Local_rho1= zeros(T,n);
    Local_rho2= zeros(T,n);
    Local_rho3=zeros(T,n);
    for ii=1:n
        Local_rho1(:,ii)= local_rho(:,1+(ii-1)*3);
        Local_rho2(:,ii)= local_rho(:,1+(ii-1)*3+1);
        Local_rho3(:,ii)= local_rho(:,1+(ii-1)*3+2);
        CC=[sum(Local_rho1,2)/n sum(Local_rho2,2)/n sum(Local_rho3,2)/n];
    end
    mm=CC(:,2);
    CV(i1,1)=std(mm)./mean(mm);
    %% Calculating the proportional abundance of each species and C.V.
end
save('ran_unsha.mat')   % Storing the results
toc