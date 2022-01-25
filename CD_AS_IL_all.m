function [COM, Nc] = CD_AS_IL(MODEL, thisnetwork,  outfile,  IsInitial, IsRefiningPartition  )
% Identifying multi-scale communities in networks by asymptotic surprise
% This is developed based on the code of Antoine Scherrer. 
% Thanks for the code of Antoine Scherrer. 
% 
% Please cite the references if you use this code: 
% [1] Vincent D. Blondel, Jean-Loup Guillaume, Renaud Lambiotte, Etienne Lefebvre
% "Fast unfolding of community hierarchies in large networks" http://arxiv.org/abs/0803.0476  
% [2] Ju Xiang, Yan Zhang, Jian-Ming Li, Hui-Jia Li, Min Li, 
% Identifying multi-scale communities in networks by asymptotic surprise, J. Stat. Mech., 2019, 2019: 033403.
% [3] Ju Xiang, Yan-Ni Tang, Yuan-Yuan Gao, Lang Liu, Yi Hao, Jian-Ming Li, Yan Zhang, Shi Chen, 
% Phase transition of Surprise optimization in community detection, Physica A, 2018, 491: 693-707.
% 
% Ju Xiang
% E-mail:xiang.ju@foxmail.com 
% 
    if nargin==0 
        MODEL=[], thisnetwork=[],  outfile=[], IsInitial=[], IsRefiningPartition=[]
        warning('!!!!!!!!!!Only for test~!!!!!!!!!!!!!!!!!!!!')
    end
    %  
    if isempty( MODEL )
        MODEL.Qualityfunction   = 'QuasiSurprise';  
        MODEL.IsWeighted        = 1; 
        MODEL.Self              = 0; 
        MODEL.Resolution        = 0.0001; 
        MODEL.MaxNiter          = [];  
        MODEL.IsDirected        = 0;
        MODEL.seed              = 111;   
    end
    %  
    if isempty( thisnetwork )
        thisnetwork.infile              = 'testnet128.wpairs';  
        thisnetwork.initialgroupfile    = ['0'];  
        thisnetwork.numnode             = -1; 
    end
    % 
    if isempty(outfile)
        outfile                         = '_detected.groups'; 
    end
    % 
    if isempty( IsInitial )
        IsInitial           = 1 ;
    end
    % 
    if isempty( IsRefiningPartition )
        IsRefiningPartition = 1; 
    end
    % % % % % % % % % % % % % % % % % % % %  
    if  IsInitial 
        thisnetwork.initialgroupfile = '1'; 
    end 
    %  
    [COM, Nc ]=LouvainFastX(MODEL, thisnetwork, outfile );            
    %  
    if  IsRefiningPartition  
        thisnetwork.initialgroupfile= '..\_detected.groups' ;MODEL.seed = [];    
        [COM, Nc ]=LouvainFastX(MODEL, thisnetwork, outfile );     
    end

end
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is developed based on code of Antoine Scherrer. 
% Thanks for the original code of Antoine Scherrer. 
% 
% Please cite these references if you use this code: 
% [1] Vincent D. Blondel, Jean-Loup Guillaume, Renaud Lambiotte, Etienne Lefebvre
% "Fast unfolding of community hierarchies in large networks" http://arxiv.org/abs/0803.0476  
% [2] Ju Xiang, Yan Zhang, Jian-Ming Li, Hui-Jia Li, Min Li, 
% Identifying multi-scale communities in networks by asymptotic surprise, J. Stat. Mech., 2019, 2019: 033403.
% [3] Ju Xiang, Yan-Ni Tang, Yuan-Yuan Gao, Lang Liu, Yi Hao, Jian-Ming Li, Yan Zhang, Shi Chen, 
% Phase transition of Surprise optimization in community detection, Physica A, 2018, 491: 693-707.
% 
% Output:
% COM: a verctor of community labels
% Q: invalid   
% Nc: the number of communities
% QUALITY: value of quality function
%  
% This is developed based on the code of Antoine Scherrer. 
% Thanks for the original code of Antoine Scherrer. 
% 
% Ju Xiang
% E-mail:xiang.ju@foxmail.com 
% 
function [COM,Nc, QUALITY] = LouvainFastX(MODEL, thisnetwork, outfile )
    ending = 0;
    if(nargin<1)   
        MODEL.Qualityfunction   = 'QuasiSurprise';  %% weighted for default, unweighted if MODEL.IsWeighted=0
        MODEL.Resolution        = 1;   % Resolution parameter  
        MODEL.IsWeighted        = [];  % 1 weighted, 0 unweighted, for  QuasiSurprise     
        MODEL.Self              = [];  % 1 weighted, 0 unweighted, for  QuasiSurprise     
        MODEL.MaxNiter          = [];  
        MODEL.recursive         = 1; 
    end
    if(nargin<2); thisnetwork.infile='netForLouvain.wpairs'; thisnetwork.initialgroupfile=[];   thisnetwork.numnode=-1;     end
    if(nargin<3); outfile='_detected.groups';end 

    if ~isfield(MODEL,'recursive') || isempty( MODEL.recursive  );  MODEL.recursive=1 ;  end   
    if isempty( MODEL.IsWeighted );   MODEL.IsWeighted    = 0;  end  %%  
    if isempty( MODEL.Self     );     MODEL.Self    = 0;  end  %%   
    if isempty( MODEL.IsDirected  );  MODEL.IsDirected    = 0;  end  %%   
    if ~isfield(MODEL,'seed') || isempty( MODEL.seed  );  MODEL.seed=123456 ;  end   
    %   

    if isfield(thisnetwork,'matrix') && ~isempty(thisnetwork.matrix)  
        M = thisnetwork.matrix;
        N = size(M ,1);     
    else
        if  ~isfield(thisnetwork,'infile') || isempty(thisnetwork.infile); error('Input is invalid!!!!!'); end  
        net0=load(thisnetwork.infile);   
        SZ=size(net0); if SZ(2)<3; net0(:,3)=1;end  %if no edge weight, let the weight of edges to be the value of 1 
        maxid_node=max(max(   net0(:,1:2) ));  
        M = sparse(net0(:,1),net0(:,2),net0(:,3),maxid_node,maxid_node); 
        N=maxid_node;  %%%%%%%%%%%%%%%
    end

    if ~MODEL.IsDirected  ; M = (M + M') ; end    %%  MODEL.IsDirected==0  

    % %
    if ~MODEL.IsWeighted  ; M= double( M>0 ) ;  end 

    %% set initial partition     
    if isfield(thisnetwork, 'initialgroup') && ~isempty(thisnetwork.initialgroup )    
        initialCOM = thisnetwork.initialgroup  ; 
        % error('thisnetwork.initialgroup is empty.');
    elseif isfield(thisnetwork, 'initialgroupfile') && ~isempty(thisnetwork.initialgroupfile )   
        switch thisnetwork.initialgroupfile
            case '0'; 	initialCOM =  [1:N];    %%%  
            case '1';   [initialCOM, COMSIZE] = Initial_partition(M,1) ;    
            otherwise % %  
                 initialGroup=load(thisnetwork.initialgroupfile);  
                 initialCOM = initialGroup(:,2);  % %  
        end
    else     
         initialCOM =  [1:N];   %%%  
    end

    if isempty( MODEL.MaxNiter );   MODEL.MaxNiter= N;  end %%  默认的迭代次序  

    global  Ln  SumLn   Nmax_for_Surprise; Ln=[];  SumLn=[];   Nmax_for_Surprise=[];     
    switch MODEL.Qualityfunction      %    calculate initial Quality of the initial partition    Et_mean        
        case 'Surprise'    %%%  
            M= M>0 ;  MODEL.Self=0; M((N+1).*[0:N-1]+1) = 0 ;  
            Nmax_for_Surprise =  200000;  
            if  N(1) <=Nmax_for_Surprise ;  N_E_max = N*(N +1) /2;
    %             Nmax_for_Surprise
                if isempty( Ln  ) || length( Ln )<N_E_max
                    Ln=log(1:N_E_max);     %%%  Ln stores the log values of all numbers from 1 to M, save the time of calculating log function. 
                    % % % SumLn=Ln;        %%$ SumLn(i) is the sum of the log values of all numbers from 1 to i, i.e., it is the log of the factorial of i  ;  
                    SumLn=zeros(1,N_E_max+1);  %% SumLn(1) SumLn(2) SumLn(3) ...SumLn(M) store the log values of   0!  1!  2! .... M! , i.e.    SumLn(i +1 ) is the log of the factorial of i  ;  
                    for i=2:N_E_max; SumLn(i   +1)=SumLn(i-1  +1)+Ln( i ) ;end 
                end
            end

        case 'QuasiSurprise'
            maxNedgeOrig =(N + MODEL.Self)*(N-1 +MODEL.Self)/2; 
            q_mean =    ( MODEL.Resolution*( N  )/(maxNedgeOrig+ eps )  )   ;
            if q_mean>1
                ending = 1; 
            end
    end  

    SZ=size(M); 
    if ~ending
        if ~MODEL.Self        
            M((N+1).*[0:N-1]+1) = 0 ;  %%消除对角
        else
            M((N+1).*[0:N-1]+1) = MODEL.Self ; 
        end   

        Q = -1;
        NODESIZE=ones(1,SZ(1));  
        NETWORK.Mnew=M;
        NETWORK.NODESIZE=NODESIZE;
        NETWORK.NnodeOrig    = N;
        NETWORK.NedgeOrig    = sum( sum(triu(M>0, ~MODEL.Self ) ) ); 
        NETWORK.maxNedgeOrig =(N + MODEL.Self)*(N-1 +MODEL.Self)/2; 
        NETWORK.LEVEL=1;
        NETWORK.Mfull=M;    
        NETWORK.Kall =  sum(M(:))  ;      %% 网络总度值
        NETWORK.km   =  NETWORK.Kall/N ;  %% 节点平均度
        NETWORK.SumWeightedEdgeOrig =NETWORK.Kall/2;    %% 所有边权之和，为总度值的一半
        NETWORK.EdgeDensity =NETWORK.NedgeOrig/(NETWORK.maxNedgeOrig+eps);    %%网络的边密度
        NETWORK.COMfull= initialCOM ;     %%% [1:SZ(1)];  
        NETWORK.initialCOM = initialCOM;  
        NETWORK.meanEdgeWeight=sum(M(:))/(sum(M(:)>0)+ 0.000000000001 ); 
        % NETWORK 
   
        [COMTY, ending] = Cluster_XJ(NETWORK, MODEL, MODEL.recursive ,  MODEL.Self) ;   
        % COMTY.COM{end} ;
        % COMTY.QUALITY(end) ; 
    % %      COM=[1:length(COMTY.COM{end});COMTY.COM{end}]';    
    else
        COMTY.COM{1} = [1:N]; 
        COMTY.MOD(1) = -1 ;  
        COMTY.QUALITY(1) = -1 ;  
    end
    Nnode = SZ(1); Nedge = sum( sum(triu(M>0, ~MODEL.Self ) ) ); output_gamma = MODEL.Resolution;     
    COM = COMTY.COM{end}; Q = COMTY.MOD(end); 
    Nc = length( unique( COMTY.COM{end} ) ) ; 
    QUALITY = COMTY.QUALITY(end) ; 
 
    fgroup_id=fopen(['..\',outfile],'w+'); fprintf(fgroup_id,'%.0f       %.0f  \r\n', [1:length(COMTY.COM{end});COMTY.COM{end}]   ); fclose(fgroup_id);   
    level    = length(COMTY.COM)-1;   
    fQmax_id=fopen('..\_fStep.Qmax','w+'); fprintf(fQmax_id,'%.0f       %f       %f      %f  \r\n',level, Q, Nc, QUALITY );  
    fclose(fQmax_id);      
     
    return ;    

end 
%%%%%%%%
function [initialCOM, COMSIZE]=Initial_partition(M, t )
    if isempty(t); t =1 ; end
    SZ= size(M); N=SZ(1); 
    initialCOM =  [1:N];  
    Ms=M*M ;   
    Ms((N+1).*[0:N-1]+1) = 0 ;  
    [i,j,v] =find(triu(M,1)); %  
    index = sub2ind(size(M),i,j); 
    v= (Ms( index )+1);
    Ms=zeros(size(M));    
	Ms(index)=v ;      % 
    Ms = Ms+Ms';   
    for j =1:t
        nodeID =randperm(N); 
    %     [s, ID] = max(CC); 
        for i =nodeID( 1: ceil(N) )
            [s, ID] = max(Ms(i,:));  
            initialCOM(i) = initialCOM(ID); 
        end
    end
    [initialCOM, COMSIZE] = reindex_com(initialCOM);
end 
 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function [COMTY, ending] = Cluster_XJ(NETWORK, MODEL, s ,self,debug,verbose)
% MODEL.Qualityfunction=  'Surprise'; 
% MODEL.Resolution=1; % Resolution parameter  
% NETWORK.Mnew      Adjecent matrix of (aggregated) network
% NETWORK.Mfull     Adjecent matrix of original full network
% NETWORK.COMfull   Community assignment of original full network    
if ~isempty( MODEL.seed  );  rng(MODEL.seed ,'twister'); MODEL.seed =MODEL.seed +111;  end   
% randperm(11) 
i=0 ;
i=i+1;if nargin < i;  error('not enough argument');end ; 
i=i+1;if nargin < i;  error('not enough argument');end ; 
i=i+1;if nargin < i;  s = 1;end ;  
i=i+1;if nargin < i;  self = 1;end ;  
i=i+1;if nargin < i;  debug = 0;end;  
i=i+1;if nargin < i;  verbose = 0;end

% % size(NETWORK.Mfull)
% % size(NETWORK.COMfull)
M = NETWORK.Mnew; 
NODESIZE = NETWORK.NODESIZE;
LEVEL = NETWORK.LEVEL;
S = size(M);
N = S(1);
Qualityfunction = MODEL.Qualityfunction; 
Resolution      = MODEL.Resolution;
MaxNiter        = MODEL.MaxNiter; 
% eps             = 0.0000000001;
ddebug = 0;
ending = 0;
eps= 10^-100; 
% Symetrize matrix taking the sum of weights
% % M = (M + M') /2;     %% 子程序外，已经对角化   
% % if (self == 0)       %% 子程序外，已经做处理
% %   M((N+1).*[0:N-1]+1) = 0;  
% % end
M2 = M;
M2((N+1).*[0:N-1]+1) = 0;   %% 

m = sum(sum(M));  
Niter = 1;

if m==0 || N == 1
    fprintf('No more possible decomposition\n');
    ending = 1;
    COMTY = 0;
    return;
end

%% Main loop
COMSTATE.Ki = sum(M); % Sum of wieght incident to super node i
COMSTATE.Ki_in = diag(M); % inner degree of super node  
% SumTot = sum(M);                                         %%%%%%%%%%%%%%%%%%%%
% SumIn = diag(M); % Sum of weight inside community i      %%%%%%%%%%%%%%%%%%%
COMSTATE.NnodeOrig    = NETWORK.NnodeOrig; %sum(NODESIZE);
COMSTATE.NedgeOrig    = NETWORK.NedgeOrig; %sum(M2(:)>0)/2;               %%%  the max number of existing edges in original network  
COMSTATE.maxNedgeOrig = NETWORK.maxNedgeOrig; %COMSTATE.NnodeOrig*(COMSTATE.NnodeOrig-1)/2;    %%%   the max number of edges in original network  
COMSTATE.EdgeDensity  =NETWORK.EdgeDensity;  %% 整个网络的边密度 
COMSTATE.SumWeightedEdgeOrig    = NETWORK.SumWeightedEdgeOrig ;  
COMSTATE.kmOrig       = NETWORK.km ;  % 
COMSTATE.Kall         = NETWORK.Kall; 
%% %% 计算每个社团的总度值，内度，节点数目
% % COM = 1:S(1); % Community of node i     set initial community assignment
COM = NETWORK.initialCOM; %   初始划分根据输入参数来确定
COMu = unique(COM);  %% 提取各个社团的标记列表
Nc =length(COMu);    %%% 社团数目
% % W_node2group  =size(N,Nc);    %%% 存储  节点 与每个社团之间的链接权重
for j=1:Nc
    Cj = find(COM==COMu(j));  %% 记录 社团j中的点
    COM( Cj ) = j ; 
    COMSTATE.KGin(j)  = full( sum(sum(M(Cj,Cj))) ) ; %% 社团内度,包含super node's inner degree （在矩阵对角线上，元素值是边权的两倍） 
    COMSTATE.KG_all(j)= sum(sum(M(Cj,:)))  ;  %% 社团总度 ,包含super node's inner degree （在矩阵对角线上，元素值是边权的两倍） 
    COMSTATE.NG(j)    =  sum( NODESIZE(Cj)  );%% 社团内节点数 （结合NODESIZE计算 原始的 节点数）   
    COMSTATE.nG(j)    =  length( Cj  );   %% 社团内节点数 （针对当前网络，不考虑NODESIZE），可能会用到
    % %     W_node2group(:,j)= sum(M(:,Cj), 2);   %% 存储 每个节点与 社团 j 中的节点Cj 之间的 链接权重总和
    % %     W_node2group(Cj,j)= W_node2group(Cj,j) - COMSTATE.Ki_in( Cj ) ;   %% 去除 节点 内度      
end
COMSTATE.M_in=   (   sum( ( COMSTATE.NG +MODEL.Self ).*(COMSTATE.NG-1  +MODEL.Self)  ) /2    );         %%%  the max number of edges within all communities of networks                    %%
COMSTATE.m_in=   (   sum(COMSTATE.KGin) /2    ) ;                 %%%  the number of existing edges within all communities of networks                   

% % sssM0 = sum(sum( M ) )
% % sssKGin = sum(sum( COMSTATE.KGin ) )
% % sssKG_all = sum(sum( COMSTATE.KG_all ) )
% % sssNG     = sum( COMSTATE.NG ) 
% % sssNG2     = sum( COMSTATE.NG(1:Nc) ) 
% % sssNODESIZE     = sum(   NODESIZE  ) 
% % ttttt=[ full(COMSTATE.KG_all)' ,  COMSTATE.KGin',   (COMSTATE.NG.*COMSTATE.NG)',  (COMSTATE.NG)'   ] 
% % SSZZ = size( ttttt )
% Ks0_in =KGin ;   Ks0_all =Ks_all; Ns0=NG;   ns0= nG;    %% 保存原始划分的统计结果，备用
% COMSTATE0=COMSTATE;  %% 保存原始划分的统计结果，备用

for k=1:N
  Neighbor{k} = find(M2(k,:));
end
%% 计算初始的  配置函数值
switch Qualityfunction      %    calculate initial Quality of the initial partition    Et_mean 
    case 'Surprise'    %  全部 不考虑自环，和权重，因为自环和权重方案还未设计好 [Surprise,Nc] = CalculateSurprise(M, M_in, m, m_in)  
        COMSTATE.M_in=   (   sum(  COMSTATE.NG .*(COMSTATE.NG-1)  ) /2    );         %%%  the max number of edges within all communities of networks                    %%
        COMSTATE.m_in=   (   sum(COMSTATE.KGin) /2    ) ;                 %%%  the number of existing edges within all communities of networks           
        [Quality ] =  CalculateSurprise(COMSTATE.maxNedgeOrig, COMSTATE.M_in,   COMSTATE.NedgeOrig    , COMSTATE.m_in);   
        
    case 'QuasiSurprise'      %[Quality,Nc ] =  CalculateQuasiSurprise(NETWORK.COMfull, NETWORK.Mfull,MODEL.IsWeighted, self);   
        COMSTATE.M_in=   (   sum( ( COMSTATE.NG +MODEL.Self ).*(COMSTATE.NG-1  +MODEL.Self)  ) /2    );         %%%  the max number of edges within all communities of networks                    %%
        COMSTATE.m_in=   (   sum(COMSTATE.KGin) /2    ) ;                 %%%  the number of existing edges within all communities of networks                   
        q_mean =   min( Resolution *(COMSTATE.M_in  )/(COMSTATE.maxNedgeOrig+ eps ) ,1)         ;     if  q_mean <=0; q_mean = 0 + eps ; end  
        q      =   (COMSTATE.m_in  + eps)/(COMSTATE.SumWeightedEdgeOrig + eps  ) ;        if  q      <=0; q      = 0  + eps ;   end  
        Quality = -COMSTATE.SumWeightedEdgeOrig*( q*log( (q_mean+ eps )/(q + eps)      + eps)   + (1-q)*log(    (1-q_mean + eps)/(1-q + eps)    + eps)    )   ;   %%直接计算 quasi-surprise
        % % 若 q_mean > q； 必须使得 q尽快的大于 q_mean  ；因为若 初始划分的 q_mean > q ，会使得 其差距愈来愈大，最终所有的社团都被拆散、分开，变成分离的点
        %         if q_mean > q;  Quality = -Quality; end     %%%   以便让程序 往 q 增大的方向行进  
        % Quality
        % Quality=0;
        % pause  
    otherwise; fprintf('no definition of the Quality function\n');
end    
 
gain = 1; 
% % % for i=1:N
% % %     indexFull_SuperNode{i}=NETWORK.COMfull==i;   %% record nodes contained in supernode i 
% % % end
 
while (gain == 1)   
    gain = 0;  
    kkk=0; 
    for i=randperm(N)   %%% 移动序列随机化
        Ci = COM(i);   kkk=kkk+1;
        %     indexCi=NETWORK.COMfull==Ci;   %% mark all nodes of Ci in original network
        NB = Neighbor{i};                  %%    AAAA
        % %     G = zeros(1,N); % Gain vector
        G = -inf(1,N); % Gain vector
        best_increase =0;  
        Quality_t=Quality; 
        QualityBest=Quality; 

        %  将点i从Ci社团移出，引起原社团改变
        COM(i) = -Ci;   
        SuperNode_Ci = COM==Ci;            % 标记属于Ci的节点        %%     AAAA
        KGin_Ci_old  = COMSTATE.KGin(Ci) ;   %% 移出前 Ci社团的内度
        KG_all_Ci_old =  COMSTATE.KG_all(Ci);  %% 移出前 Ci社团的总度
        NG_Ci_old = COMSTATE.NG(Ci) ;          %% 移出前 Ci社团的节点总数
        nG_Ci_old = COMSTATE.nG(Ci) ;          %% 移出前 Ci社团的super节点总数 
        % 更新原Ci社团信息
        COMSTATE.KGin(Ci)  = COMSTATE.KGin(Ci) -  (2*sum(M(SuperNode_Ci,i)) +    M(i,i) );  %% 新的内度， 内度减小
        
        COMSTATE.KG_all(Ci)= COMSTATE.KG_all(Ci) - COMSTATE.Ki(i);             %% 社团总度减小      
        COMSTATE.NG(Ci)    =  COMSTATE.NG(Ci)-  NODESIZE(i);      %% 社团内节点数减少  
        COMSTATE.nG(Ci)    =  COMSTATE.nG(Ci)- 1 ;      %% 社团内节点数减少 

        % %      W_node2group(:,Ci)=  W_node2group(:,Ci) -M(:,i);   %% 存储 每个节点与 社团 j 中的节点Cj 之间的 链接权重总和

        Cnew = Ci;    Cnew_t=Cnew; 
        %     COMnew=COM;   %%%%%    
        %     COM(i) = -1;                         
        % %     CNj1 = find(COM==Ci);
        % %     SumTot(Ci) = SumTot(Ci) - K(i);
        % %     SumIn(Ci) = SumIn(Ci) - 2*sum(M(i,CNj1)) - M(i,i);
            %%%%%%%%%%%%%%%%%
        %     tic
        for j=1:length(NB)
        
            Cj = COM(NB(j));
            % %       if (G(Cj) == 0)     %% G(Cj) == 0  表明向该社团的移动还没有尝试过 
            if (G(Cj) == -inf)     %% G(Cj) == -inf  表明向该社团的移动还没有尝试过 
                % %         CNj = find(COM==Cj);        
                % %         Ki_in = 2*sum(M(i,CNj));
                % %         G(Cj) = Ki_in/m - 2*K(i)*SumTot(Cj)/(m*m);  

                %  将点i放入Cj社团，引起Cj社团改变 
                SuperNode_Cj  = COM==Cj;            % 标记属于Cj的节点 
                KGin_Cj_new   = COMSTATE.KGin(Cj) +  2*sum(M(SuperNode_Cj,i)) +M(i,i);     %%Cj社团的新内度， 内度增加
                KG_all_Cj_new = COMSTATE.KG_all(Cj) + COMSTATE.Ki(i);             %% 社团Cj的总度增加     
                NG_Cj_new = COMSTATE.NG(Cj) + NODESIZE(i) ;     %% 社团内节点数增加   
                 
                Is_Pg_LT_P = 0;
                switch Qualityfunction    %%  calculate Quality after moving node i into community Cj        
                    case 'Surprise' ;
                        %                 M_in=   ceil(   sum(  NG_t.*(NG_t-1)  ) /2    );        %%%  the max number of edges within all communities of networks                    %%
                        %                 m_in=   ceil(   sum(KGin_t) /2    );                 %%%  the number of existing edges within all communities of networks                

                        M_in_t=  (  COMSTATE.M_in -  NG_Ci_old*(NG_Ci_old-1)/2  -  COMSTATE.NG(Cj)*(COMSTATE.NG(Cj)-1)/2   +  COMSTATE.NG(Ci)*(COMSTATE.NG(Ci)-1)/2   +  NG_Cj_new*(NG_Cj_new-1)/2   );        %%%  the max number of edges within all communities of networks                    %%
                        m_in_t=  (  COMSTATE.m_in -  KGin_Ci_old/2   -COMSTATE.KGin(Cj)/2    + COMSTATE.KGin(Ci)/2  +  KGin_Cj_new/2     );                 %%%  the number of existing edges within all communities of networks           
                        %                 COMSTATE.M_in_t= M_in_t;
                        %                 COMSTATE.m_in_t = m_in_t

                        [Quality_t ] =  CalculateSurprise(COMSTATE.maxNedgeOrig, M_in_t,   COMSTATE.NedgeOrig    , m_in_t);    
                        %  toc
                    case 'QuasiSurprise' ;
                        M_in_t=  (  COMSTATE.M_in -  (NG_Ci_old +MODEL.Self)*(NG_Ci_old-1 +MODEL.Self)/2  -  (COMSTATE.NG(Cj) +MODEL.Self)*(COMSTATE.NG(Cj)-1 +MODEL.Self)/2   +  (COMSTATE.NG(Ci) +MODEL.Self)*(COMSTATE.NG(Ci)-1 +MODEL.Self)/2   +  (NG_Cj_new +MODEL.Self)*(NG_Cj_new-1 +MODEL.Self)/2   );        %%%  the max number of edges within all communities of networks                    %%
                        m_in_t=  (  COMSTATE.m_in -  KGin_Ci_old/2   -COMSTATE.KGin(Cj)/2    + COMSTATE.KGin(Ci)/2  +  KGin_Cj_new/2     );                 %%%  the number of existing edges within all communities of networks           
                        % %             q_mean =   ( M_in_t   )/(COMSTATE.maxNedgeOrig+ eps  )     ;   q_mean =q_mean^Resolution;         if  q_mean <=0; q_mean = 0 + eps ; end  
                        q      =   (m_in_t  + eps )/(COMSTATE.SumWeightedEdgeOrig + eps ) ;        if  q      <=0; q      = 0  + eps ;   end  
                        q_mean =  Resolution * ( M_in_t   )/(COMSTATE.maxNedgeOrig+ eps  );
                        if q <=q_mean  && q_mean>=1  
                            q0 = (COMSTATE.m_in  + eps )/(COMSTATE.SumWeightedEdgeOrig + eps ) ; 
                            if q < q0 
                                Is_Pg_LT_P = 1 ;   % 针对(多分辨率方法的)极端情况 
                            else
                                Is_Pg_LT_P = 0;
                            end
                        else
                            Is_Pg_LT_P = 0;
                        end
                        %         
                        if Is_Pg_LT_P == 0 
                            q_mean =   max(eps, min( q_mean ,1)  )   ;           
                            if  q_mean <=0; q_mean = 0 + eps ; end  
                            %                 q_mean =   (M_in_t )/(COMSTATE.maxNedgeOrig+ eps )     ;            if  q_mean <=0; q_mean = 0 + eps ; end  
                            %                 q      =   (m_in_t  )/(COMSTATE.SumWeightedEdgeOrig + eps ) ;        if  q      <=0; q      = 0  + eps ;   end  
                            Quality_t = -COMSTATE.SumWeightedEdgeOrig*( q*log(  (q_mean + eps)/(q + eps)      + eps)   + (1-q)*log(    (1-q_mean + eps)/(1-q + eps)    + eps)    )   ;   %%直接计算 quasi-surprise
                            % % 若 q_mean > q； 必须使得 q尽快的大于 q_mean  ；因为若 初始划分的 q_mean > q ，会使得 其差距愈来愈大，最终所有的社团都被拆散、分开，变成分离的点
                            %                 if q_mean > q;  Quality_t = -Quality_t; end     %%%   以便让程序 往 q 增大的方向行进              
                            %         COMSTATE
                            %         [M_in_t    m_in_t    q_mean*1000     q*1000    Quality_t    Resolution ]
                            %         if   q_mean <0  ||  q<0; fprintf('q_mean*1000=%f    q*1000=%f \n', q_mean*1000 ,     q*1000  );   pause ; end
                        end 
                end        
                
                if Is_Pg_LT_P
                    G(Cj) = -10^100 ;   % 不允许 移动 
                else                    
                    G(Cj) = Quality_t- Quality;   %% increase of Quality
                end
                
                %         Quality_t  
                if (ddebug)
                    fprintf('Gain for comm %d => %g\n',Cj,G(Cj));
                end
                if G(Cj) > best_increase 
                    best_increase = G(Cj);
                    QualityBest=Quality_t ;
                    Cnew_t = Cj;
                end
            end
        end
     
        % %     pause; 
        if best_increase > 0
            Cnew = Cnew_t;
            if (debug)
                fprintf('Move %d => %d\n',i-1,Cnew-1);
            end
            %       Cost(i) = best_increase;
        end
        %     Ck = find(COM==Cnew);
        %     SumIn(Cnew) = SumIn(Cnew) + 2*sum(M(i,Ck));
        %     SumTot(Cnew) = SumTot(Cnew) + K(i);
            %  将点i最终放入（产生最大增量的）Cnew 社团，引起Cj=Cnew社团改变(最终的更新) 
        Cj=Cnew;
        SuperNode_Cj = COM==Cj;            % 标记属于Ci的节点 
        KGin_Cj_new  = COMSTATE.KGin(Cj)  +  (2*sum(M(SuperNode_Cj,i)) +    M(i,i) ); %% 新的内度， 内内度增加
        KG_all_Cj_new = COMSTATE.KG_all(Cj) + COMSTATE.Ki(i);             %% 社团Cj的总度增加         
        NG_Cj_new = COMSTATE.NG(Cj) +  NODESIZE(i);     %% 社团内节点数增加   
        nG_Cj_new = COMSTATE.nG(Cj) +  1;     %% 社团内节点数增加   

        COMSTATE.M_in=  (  COMSTATE.M_in -  (NG_Ci_old +MODEL.Self)*(NG_Ci_old-1 +MODEL.Self)/2  -  (COMSTATE.NG(Cj) +MODEL.Self)*(COMSTATE.NG(Cj)-1 +MODEL.Self)/2   +  (COMSTATE.NG(Ci) +MODEL.Self)*(COMSTATE.NG(Ci)-1 +MODEL.Self)/2   +  (NG_Cj_new +MODEL.Self)*(NG_Cj_new-1 +MODEL.Self)/2   );        %%%  the max number of edges within all communities of networks                    %%
        COMSTATE.m_in=  (  COMSTATE.m_in -  KGin_Ci_old/2   -COMSTATE.KGin(Cj)/2    + COMSTATE.KGin(Ci)/2  +  KGin_Cj_new/2     );                    %%%  the number of existing edges within all communities of networks           

        %更新新社团的信息
        COMSTATE.KGin(Cj)  = KGin_Cj_new;         %% 必须更新
        COMSTATE.KG_all(Cj)= COMSTATE.KG_all(Cj) + COMSTATE.Ki(i);            %% 社团总度增加
        COMSTATE.NG(Cj)    =  NG_Cj_new   ;   %% 必须更新
        COMSTATE.nG(Cj)    =  nG_Cj_new;                    %% 社团内节点数增加         
   
        COM(i) = Cnew;   %% 最终更新  COM   列表  
        %     NETWORK.COMfull(indexFull_SuperNode{i}) = Cnew;  
        Quality=QualityBest ;    %  Quality 
        %     if ~isreal(Quality);
        %       COMSTATE
        %       pause
        %     end
        %     ttttt=sum(  COMSTATE.NG.*(COMSTATE.NG-1)/2) 
        %     wwww=COMSTATE.M_in
        %     COMSTATE 
        %     pause 
        %     q.COM1=COM;
        %     q.COM2=NETWORK.COMfull;
        %     q
        %     pause 
        %     Quality=Quality_t ;    % update community assignment and Quality
        if (Cnew ~= Ci) ;    gain = 1;    end

        % kk= sum( COMSTATE.KG_all  )   
        % nn=sum( COMSTATE.NG  )  
        % pause
    end
    %   toc

 
  %%
    %   sCost = sum(Cost);
    if (debug)
        [C2 S2] = reindex_com(COM);
        Nco = length(unique(COM));
        Nco2 = length(S2(S2>1));
        mod = CalculateModularityCM(COM,M,Resolution);
        fprintf('It %d - Mod=%f %d com (%d non isolated)\n',Niter,mod,Nco,Nco2);
    end
  
	Niter = Niter + 1 ; 
    %     q.COM1=COM;
    %     q.COM2=NETWORK.COMfull;
    %     q
    if MaxNiter~=-1 && Niter > MaxNiter; break; end     
  
end

Niter = Niter - 1;
[COM, COMSIZE] = reindex_com(COM);
COMTY.COM{1} = COM;
COMTY.SIZE{1} = COMSIZE;
COMTY.MOD(1) =  -1;
 
COMTY.QUALITY(1) = full( Quality ) ;  %  CalculateModularity(COM,M);   %% 避免重复计算 Quality
COMTY.Niter(1) = Niter;
% COMTY
% pause
switch Qualityfunction    %% for Qualityfunction of  searching for minima                
    case 'Hamijmmmmm';  COMTY.QUALITY(1) = - Quality  ;           
end
 
 
%% Perform part 2
if (s == 1)
  
    Mnew = M;
    Mold = Mnew;
    COMcur = COM;
    COMfull = COM;

    k = 2;
    if (debug)
        Nco2 = length(COMSIZE(COMSIZE>1));      fprintf('Pass number 1 - %d com (%d iterations)\n',Nco2,Niter);
    end
	while 1
        Mold = Mnew;
        S2 = size(Mold);
        Nnode = S2(1);
    %     NODESIZEnew=[];

        COMu = unique(COMcur);
        Ncom = length(COMu);
        ind_com = zeros(Ncom,Nnode);
        ind_com_full = zeros(Ncom,N);

        for p=1:Ncom
            ind = find(COMcur==p);
            ind_com(p,1:length(ind)) = ind;

            ind = find(COMfull==p);
            ind_com_full(p,1:length(ind)) = ind;
            %       NODESIZEnew(p)=length(ind);   %% number of nodes in each supernode of new network
        end
    
        Mnew = zeros(Ncom,Ncom);
        for m=1:Ncom
            for n=m:Ncom
                ind1 = ind_com(m,:);
                ind2 = ind_com(n,:);
                Mnew(m,n) = sum(sum(Mold(ind1(ind1>0),ind2(ind2>0))));
                Mnew(n,m) = Mnew(m,n);
                %         Mnew(n,m) = sum(sum(Mold(ind1(ind1>0),ind2(ind2>0))));
            end
        end
        % %         sssM0 = sum(sum(M) ) 
        % %         sssMnew = sum(sum(Mnew) )
        % %         [ sum(Mnew,2), diag(Mnew),COMSIZE(:).*COMSIZE(:),   COMSIZE(:) ] 
        
        %     NODESIZEnew
            NETWORK.Mnew=Mnew;
        %     NETWORK.Mfull=NETWORK.Mfull;  %% unchange   
            NETWORK.COMfull=COMfull;  
            NETWORK.initialCOM = [1:Ncom] ;   %%%%%%%%%%%新加入初始划分的设置%%  对于新的 supernode 网络，初始化 每个 supernode 作为一个社团 
            NETWORK.LEVEL=LEVEL+1; 
            NETWORK.NODESIZE=COMSIZE;  
        %     szzz=size(NETWORK.COMfull) 
        %     NETWORK
        %     pause
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 此处 是否必须要考虑自环 self =1 ！若不考虑自环，矩阵中的节点内度将丢失！！！！！！
        [COMt, e] = Cluster_XJ(NETWORK, MODEL, 0,1,debug,verbose);    %   no iteration, so s=0
        %     Ncom =Ncom
    
        if (e ~= 1)
            COMfull = zeros(1,N);
            COMcur = COMt.COM{1};
            for p=1:Ncom
                ind1 = ind_com_full(p,:);
                COMfull(ind1(ind1>0)) = COMcur(p);
            end
            % XXXXX [COMfull, COMSIZE] = reindex_com(COMfull);  % % XXXX问题  可能改变了 社团编号顺序         
            COMSIZE = sparse(COMcur,1, COMSIZE)  ;
            % % COMSIZE = sparse(COMfull,1, 1)  ;            
            COMTY.COM{k} = COMfull;
            COMTY.SIZE{k} = COMSIZE;
            COMTY.MOD(k) = COMt.MOD(1);   
            COMTY.QUALITY(k) = COMt.QUALITY(1)  ;    
            COMTY.Niter(k) = COMt.Niter;
            if (debug)
                Nco2 = length(COMSIZE(COMSIZE>1));
                fprintf('Pass number %d - %d com\n',k,Nco2);
            end
            %%   
% %             Ind = (COMfull == COMTY.COM{k-1});   
% %             if (sum(Ind) == length(Ind))
            if all(  COMfull == COMTY.COM{k-1}  )
                if (debug)
                    fprintf('Identical segmentation => End\n');
                end
                return;
            end
        else
            if (debug);         fprintf('Empty matrix => End\n');        end
            return;
        end
    
        k = k + 1;
	end
end

end
 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%% Re-index community IDs
function [C Ss] = reindex_com(COMold)
    C = zeros(1,length(COMold));
    COMu = unique(COMold);
    S = zeros(1,length(COMu));
    for l=1:length(COMu)
        S(l) = length(COMold(COMold==COMu(l)));
    end
    [Ss INDs] = sort(S,'descend');

    for l=1:length(COMu)
        C(COMold==COMu(INDs(l))) = l;
    end

end

%% Compute ModulartiyCM
function [MOD,Nc] = CalculateModularityCM(C,Mat,Resolution)
m = sum(sum(Mat));
MOD = 0;
COMu = unique(C); 
Nc =length(COMu);  
for j=1:length(COMu)
    Cj = find(C==COMu(j));
    Ec = sum(sum(Mat(Cj,Cj)));
    Et = sum(sum(Mat(Cj,:)));
    if Et>0
        MOD = MOD + Ec/m-Resolution*(Et/m)^2;
    end
end
end
 

  
%% Compute Surprise of a given partition, noweight
% function [Surprise,Nc] = CalculateSurprise(N, nG, KG , KGin  )  
% (COMSTATE.maxNedgeOrig, M_in_t,   COMSTATE.NedgeOrig    , m_in_t); 
function [Surprise] = CalculateSurprise(M, M_in, m, m_in)   
% no weight, no selfloop: without considering edge weigth and self-loop
% com  contains community ID of each node
% Mat  Matrix give full unweigted network
% NODESIZE record the number of nodes in each supernode 
% self= 0;    %%%    
%% translate the multiplication into being addition by log function  
global    SumLn   Nmax_for_Surprise;    %%% 如果体系较小，提前计算  SumLn 减少重复计算
% Nmax=max( 5000, Nmax_for_Surprise ); 
if      M<=Nmax_for_Surprise*(Nmax_for_Surprise+1 )/2 
    % tic 
    % global Ln    SumLn;   
    % a='ssssss'
    % Ln=log(1:M);     %%%  Ln stores the log values of all numbers from 1 to M, save the time of calculating log function. 
    % % % % SumLn=Ln;        %%$ SumLn(i) is the sum of the log values of all numbers from 1 to i, i.e., it is the log of the factorial of i  ;  
    % SumLn=zeros(1,M+1);  %% SumLn(1) SumLn(2) SumLn(3) ...SumLn(M) store the log values of   0!  1!  2! .... M! , i.e.    SumLn(i +1 ) is the log of the factorial of i  ;  
    % for i=2:M; SumLn(i   +1)=SumLn(i-1  +1)+Ln( i ) ;end 

    % SumLn
    Surprise=0;    
    logdenominator= SumLn(M    +1)- SumLn(m  +1)-SumLn(M-m  +1);  

    % for i=m_in:min(M_in,m) 
    %     Surprise= Surprise + exp(  [SumLn(M_in  +1)-SumLn(i   +1)-SumLn(M_in-i  +1)] + [ SumLn(M-M_in  +1)-SumLn(m-i  +1)-SumLn( (M-M_in) - (m-i)   +1) ] - logdenominator    );   
    %     % %  
    % end
    % i=m_in-min(M_in,m) 
    %%  
    % k=0;    
    % for i=m_in:min(M_in,m) 
    %     k=k+1;      X(k)=(SumLn(M_in  +1)-SumLn(i   +1)-SumLn(M_in-i  +1) )  +  ( SumLn(M-M_in  +1)-SumLn(m-i  +1)-SumLn( (M-M_in) - (m-i)   +1) )    - logdenominator   ;   
    % end
    % 
    i=m_in:min(M_in,m) ;   X = 0; 
% %     X=(SumLn(M_in  +1)-SumLn(i   +1)-SumLn(M_in-i  +1) )  +  ( SumLn(M-M_in  +1)-SumLn(m-i  +1)-SumLn( (M-M_in) - (m-i)   +1) )    - logdenominator   ;   
    X= SumLn(M_in  +1)-SumLn(i   +1)-SumLn(M_in-i  +1)      +   SumLn(M-M_in  +1)-SumLn(m-i  +1)-SumLn( (M-M_in) - (m-i)   +1)       - logdenominator   ;   

    % com
    % sssss
    % M
    % M_in
    % m
    % m_in
    % % disp('SaveSaveSaveSaveSaveSaveSaveSaveSave ')

    [Xmax,ii]=max(X)  ;
    Surprise= - (  log(  sum( exp( X-Xmax  ) )  ) + Xmax ) ;
    % toc 
else
    %%
% % %     tic 
    sumLn_M=0;   for j=1:M;    sumLn_M=sumLn_M+log(j);  end;
    sumLn_m=0;   for j=1:m;    sumLn_m=sumLn_m+log(j);  end;
    sumLn_M_m=0; for j=1:(M-m);  sumLn_M_m=sumLn_M_m+log(j);end; 
    %
    sumLn_M_in=0;        for j=1:M_in;    sumLn_M_in=sumLn_M_in+log(j);  end;
    sumLn_m_in=0;        for j=1:m_in;    sumLn_m_in=sumLn_m_in+log(j);  end;
    sumLn_M_in_m_in=0;   for j=1:(M_in-m_in);    sumLn_M_in_m_in=sumLn_M_in_m_in+log(j);  end;
    %
    sumLn_M_M_in=0;      for j=1:(M-M_in);    sumLn_M_M_in=sumLn_M_M_in+log(j);  end;
    sumLn_m_m_in=0;      for j=1:(m-m_in);    sumLn_m_m_in=sumLn_m_m_in+log(j);  end;
    sumLn_M_M_in_m_m_in=0;   for j=1:((M-M_in)-(m-m_in));    sumLn_M_M_in_m_m_in=sumLn_M_M_in_m_m_in+log(j);  end;
    %
    k=1; X(k)=(sumLn_M_M_in-sumLn_m_m_in-sumLn_M_M_in_m_m_in)+  ( sumLn_M_in- sumLn_m_in- sumLn_M_in_m_in) - (sumLn_M- sumLn_m- sumLn_M_m ); 
    for j=m_in+1: min(M_in,m) 
        k=k+1; X(k)= X(k-1) -log(j)+log( M_in- (j-1) )  + log(m- (j-1)) - log(  (M-M_in)-(m- (j) )  ) ; 
    end

    [Xmax,ii]=max(X) ;
    Surprise= - (  log(  sum( exp( X-Xmax  ) )  ) + Xmax ) ;
end
 
 
end

%% %% Compute Quasi-Surprise of a given partition 
function [QuasiSurprise,Nc] = CalculateQuasiSurprise(com,Mat, IsWeighted  , self, NODESIZE)   
% Quasi Surprise see PRE 92, 022816(2015)  
% Detecting communities using asymptotical surprise
% no weight, no selfloop: without considering edge weigth and self-loop
% com  contains community ID of each node
% Mat  Matrix give full unweigted network
% IsWeighted : 1 consider weight      0 without weight, IsWeighted=1 for default  
% NODESIZE record the number of nodes in each supernode 
%%%  Nc = length(unique(com));  
    Fself = ~self; 
    [com Ss] = reindex_com(com); 
    eps=0.00000000000000001; 
    maxid_group=max( com ); 
    SZ=size(Mat); maxid_node=SZ(1);
    [i,j,w]=find(   triu(  Mat, Fself  )  )  ;   %% aviod selfloop 
    if ~IsWeighted ;  w(:)=1;  end %% to be unweighted network
% %     meanEdgeWeight=(sum(w) )/(length(w)+eps ) ;    
    NODESIZE=ones(1,maxid_node);    % NODESIZE  records the real number of nodes in each supernode , here the full network is given 
    net0 =[i,j,w];
    meanEdgeWeight=1;   %%%改变此处不意味着适用于加权网络%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     (com)
%     [( net0  ), net0(:,3)  ]
    net1 =[com( net0(:,1:2) ), net0(:,3)];      %%%graph with group ID  
%     net1 =[com( i )', com(j)',  w];       %%%graph with group ID     
    netG= sparse( [net1(:,1);net1(:,2)],  [net1(:,2);net1(:,1)]  ,[net1(:,3);net1(:,3)],maxid_group, maxid_group )  ;
    KG= sum(netG) ;
    KGin= diag(netG); 
%     nG=sparse( com,1,1, maxid_group,1 );    %% nG denotes the number of nodes in each Group. 
    nG=sparse( com,1,NODESIZE, maxid_group,1 );    %% nG denotes the number of nodes in each Group.  realistic number of nodes in groups by NodeSize 
    Ksum=sum(KG);
M   =   floor( maxid_node*(  maxid_node- Fself  )/2 );      %%%  the max number of edges in network with maxid_node nodes
M_in=   ceil(   sum(  nG.*(nG-Fself)  ) /2    )  ;     %%%  the max number of edges within all communities of networks 
m   =   ceil(   Ksum/2  )     ;                    %%%  the number of existing edges in networks 
m_in=   ceil(   sum(KGin) /2    )  ;              %%%  the number of existing edges within all communities of networks 

% Ksum
% sum(KGin)
% full( netG )
% stop
    q_mean = M_in/(M+ eps)     ;    if  q_mean <=0; q_mean = q_mean + eps ; end  
    q      = m_in/(m + eps) ;       if  q      <=0; q      = q  + eps ;   end  
    QuasiSurprise = -m*( q*log( q_mean/(q + eps)      + eps)   + (1-q)*log(    (1-q_mean)/(1-q + eps)    + eps)    )         ;
    QuasiSurprise=full(QuasiSurprise);
    Nc=maxid_group;
    %% 若 q_mean > q； 必须使得 q尽快的大于 q_mean  ；因为若 启动点 q_mean > q ，会使得 其差距愈来愈大，最终所有的社团都被分开，变成分离的点
   if q_mean > q; QuasiSurprise = -QuasiSurprise; end     %%%
% % % stop
end




