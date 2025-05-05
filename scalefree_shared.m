clear
clc
tic
%% setting  parameters
n=16; % Coarse-scale patch size
s=3; % Species number
num=1024;% Fine-scale patch number
T = 1500;% generations
H = [0 1 0;0 0 1;1 0 0]; % competition matrix
i2=-10;
i3=length(i2);
CV=zeros(i3,1);
load('scal.mat') %  scaleefree network matrix
network=M;
load('scala.mat')% Fine-scale patch  matrix
for i1=1:i3
    im=i1+100;
    io= i2(i1);
    lamida = 2^(io) ; %dispersal probability D
    local_rho=zeros(T,s*n);
    load('initialspecies.mat');
    for rep = 1:T
        rep
        for t=1:num
            focous_p = zeros(n,1);    %selecting empty patch
            for ig=1:n
                focous_p(ig,:) = randperm(num,1)+(ig-1)*num;  % focal patch
            end
            I_species(focous_p)=0;
            species=I_species;
            pd1 = rand(n,1);
            pd = ones(n,1);pd(pd1<lamida)=0;
            for sc =1:n  %occupation patch
                %%
                if pd(sc)==1
                    M=a{sc,1};
                    part_p=focous_p(sc)-(sc-1)*num;
                    link=find(M(:,part_p)==1);
                    if length(link)==1
                        species(focous_p(sc)) =I_species(link,sc);
                    else
                        random_num = link(randperm(numel(link),2));
                        s1= I_species(random_num,sc);
                        s11=s1(1,1);
                        s12=s1(2,1);
                        san = H(s11,s12)<rand(1);
                        species(focous_p(sc))= s1(1+san);
                    end
                    %% selecting species within Fine-scale patch occupy empty patch
                    %%
                else
                    linju = find(network(:,sc)==1);
                    link_local_cell=I_species(:,linju);
                    link_local_cell( link_local_cell==0)=[];
                    L1=length( link_local_cell);
                    x=randperm(L1,2);
                    competitotrs_exs=link_local_cell(x);
                    HH2=H( competitotrs_exs(1), competitotrs_exs(2));
                    species(focous_p(sc))=  competitotrs_exs(2-HH2);
                    %% selecting species within Coarse-scale patch occupy empty patch
                end
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
save('scal_sha.mat')  %Store the results
toc