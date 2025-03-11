function  [pS,T,W] = channel_discrimination_5copies_primal_mine(C,protocol,varargin)
% This function implements the primal SDP for discriminating a set of N
% quantum channels {C_i}_i when k=5 copies are available.
%
% The function outputs the optimal succes probability pS,
% the set of testers T(:,:,1), ..., T(:,:,N) that attains pS,
% and the process W = \sum_i T(:,:,i)
% If causally separable testeres are considered, W is a tensor which stores the process from
% 1 to 2 and the process from 2 to 1. 
% That is, W(:,:,1)=W12, and W(:,:,2)=W21.
%
% The set of channels shoud be storered in a variable C
% For instance: C(:,:,1)=C1, C(:,:,2)=C2, ..., C(:,:,N)=CN  
%
% When protocol==1, Paralel protocols are considered
% When protocol==2, Sequential protocols are considered
% When protocol==3, Separable protocols are considered
% When protocol==4, General protocols are considered
%
% If the channels have different input/output dimension, one should write
% channel_discrimination_5copies_primal(C,protocol,[dIn dOut])
% If the distribution of the channels is not the uniform one, i.e., p_i=1/N, one should write
% channel_discrimination_2copies_primal(C,protocol,[dIn dOut], p_i)
% where p_i is an array with the distribution p_i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=size(C,3); %Obtain the number of channels N
k=5; %Set the number of uses k equals 5

% Switch that analyses the extra inputs varargin
switch length(varargin)
    case 0  %If no extra input is given, assume dIn=dOut and uniform p_i
        d=sqrt(size(C(:,:,1),1));
        dIn=d;
        dOut=d;
        DIM=[d d d d d d d d d d]; %%ADDED [... d d]
        p_i=ones(1,N)/N;
    case 1 %If one extra input is given, assume uniform p_i
        dIn=varargin{1}(1);
        dOut=varargin{1}(2);
        DIM=[dIn dOut dIn dOut dIn dOut dIn dOut dIn dOut]; %%ADDED [... dIn dOut]
        p_i=ones(1,N)/N;
    case 2 %If two extra inputs are given, use information from extra inptus
        dIn=varargin{1}(1);
        dOut=varargin{1}(2);
        DIM=[dIn dOut dIn dOut dIn dOut dIn dOut dIn dOut];  %%ADDED [... dIn dOut]
        p_i=varargin{2};
end

d=2; %%THIS ONLY APPEATS IN k=3, WHY?
% Switch that treats the 4 possible different protocols
switch protocol   
%%%%%%%%%%%%%%%%%%%%% PARALLEL TESTERS %%%%%%%%%%%%%%%%%%%%%
    case 1
cvx_begin SDP
%cvx_solver sedumi
     variable T(dIn^(2*k) dOut^(2*k) N) complex semidefinite
     expression pS ;
     expression W ;
     pS=0;
     W=0;     
     for i=1:N
         pS = pS + trace(p_i(i)*T(:,:,i) * Tensor(C(:,:,i),k) );
         W = W + T(:,:,i);
     end
     
     W == TR(W,[2 4 6 8 10],DIM);
     trace(W) == dOut^k;
  
     maximise pS;  
cvx_end
%%%%%%%%%%%%%%%%%%%%% SEQUENTIAL TESTERS %%%%%%%%%%%%%%%%%%%%%
    case 2
cvx_begin SDP
%cvx_solver sedumi
     variable T(dIn^(2*k) dOut^(2*k) N) complex semidefinite
     expression pS ;
     expression W ;
     pS=0;
     W=0;
     
     for i=1:N
         pS = pS + trace(p_i(i)*T(:,:,i) * Tensor(C(:,:,i),k) );
         W = W + T(:,:,i);
     end
     
     W == TR(W,[10],DIM); %%COND A
     PartialTrace(W,[10 9],[d d d d d d d d d d]) == kron(PartialTrace(W,[10 9 8],[d d d d d d d d d d]),eye(d)/d); %% COND B
     PartialTrace(W,[10 9 8 7],[d d d d d d d d d d]) == kron(PartialTrace(W,[10 9 8 7 6],[d d d d d d d d d d]),eye(d)/d); %% COND B
     PartialTrace(W,[10 9 8 7 6 5],[d d d d d d d d d d]) == kron(PartialTrace(W,[10 9 8 7 6 5 4 ],[d d d d d d d d d d]),eye(d)/d); %% COND B
     PartialTrace(W,[10 9 8 7 6 5 4 3],[d d d d d d d d d d]) == kron(PartialTrace(W,[10 9 8 7 6 5 4 3 2],[d d d d d d d d d d]),eye(d)/d); %% COND B
     trace(W) == dOut^k; %% MODIFIED TO d^4
  
     maximise pS;  
cvx_end 

%Conditions for k=2 and k=3
% W == TR(W,[4],DIM);
%      TR(W,[4 3],DIM) == TR(W,[4 3 2],DIM);
%      trace(W) == dOut^2;

% W == TR(W,[6],DIM);
 %     PartialTrace(W,[6 5],[d d d d d d]) == kron(PartialTrace(W,[6 5 4],[d d d d d d]),eye(d)/d);
 %     PartialTrace(W,[6 5 4 3],[d d d d d d]) == kron(PartialTrace(W,[6 5 4 3 2],[d d d d d d]),eye(d)/d);
 %     trace(W) == dOut^3;

%%%%%%%%%%%%%%%%%%%%% GENERAL TESTERS %%%%%%%%%%%%%%%%%%%%%
    case 4
cvx_begin SDP
%cvx_solver sedumi
     variable T(dIn^(2*k) dOut^(2*k) N) complex semidefinite
     expression pS ;
     expression W ;
     pS=0;
     W=0;
     
     for i=1:N
         pS = pS + trace(p_i(i)*T(:,:,i) * Tensor(C(:,:,i),k) );
         W = W + T(:,:,i);
     end
    %%% Marginals
%L_A(W)=W 
    TR(W,[1 2 3 4 5 6 7 8],DIM) == TR(W,[1 2 3 4 5 6 7 8 10],DIM); %L_E(W)=W 
    TR(W,[1 2 3 4 5 6 9 10],DIM) == TR(W,[1 2 3 4 5 6 8 9 10],DIM); %L_D(W)=W
    TR(W,[1 2 3 4 7 8 9 10],DIM) == TR(W,[1 2 3 4 6 7 8 9 10],DIM); %L_CW)=W 
    TR(W,[1 2 5 6 7 8 9 10],DIM) == TR(W,[1 2 4 5 6 7 8 9 10],DIM);%L_B(W)=W 
    TR(W,[3 4 5 6 7 8 9 10],DIM) == TR(W,[2 3 4 5 6 7 8 9 10],DIM);%L_A(W)=W 
    


%L_AB(W)=W 
    TR(W,[5 6 7 8 9 10],DIM) + TR(W,[2 4 5 6 7 8 9 10],DIM) == TR(W,[2 5 6 7 8 9 10],DIM) + TR(W,[4 5 6 7 8 9 10],DIM); %L_AB(W)=W
    TR(W,[1 2 3 4 9 10],DIM) + TR(W,[1 2 3 4 6 8 9 10],DIM) == TR(W,[1 2 3 4 6 9 10],DIM) + TR(W,[1 2 3 4 8 9 10],DIM); %L_CD(W)=W 
    TR(W,[1 2 5 6 9 10],DIM) + TR(W,[1 2 4 5 6 8 9 10],DIM) == TR(W,[1 2 4 5 6 9 10],DIM) +  TR(W,[1 2 5 6 8 9 10],DIM); %L_BD(W)=W 
    TR(W,[1 2 7 8 9 10],DIM) + TR(W,[1 2 4 6 7 8 9 10],DIM) == TR(W,[1 2 4 7 8 9 10],DIM) +  TR(W,[1 2 6 7 8 9 10],DIM); %%L_BC(W)=W 
    TR(W,[3 4 5 6 9 10],DIM) + TR(W,[2 3 4 5 6 8 9 10],DIM) == TR(W,[2 3 4 5 6 9 10],DIM) +  TR(W,[3 4 5 6 8 9 10],DIM); %L_AD(W)=W 
    TR(W,[3 4 7 8 9 10],DIM) + TR(W,[2 3 4 6 7 8 9 10],DIM) == TR(W,[2 3 4 7 8 9 10],DIM) +  TR(W,[3 4 6 7 8 9 10],DIM); %L_AC(W)=W 
    TR(W,[3 4 5 6 7 8],DIM) + TR(W,[2 3 4 5 6 7 8 10],DIM) == TR(W,[2 3 4 5 6 7 8],DIM) +  TR(W,[3 4 5 6 7 8 10],DIM); %L_AE(W)=W
    TR(W,[1 2 5 6 7 8],DIM) + TR(W,[1 2 4 5 6 7 8 10],DIM) == TR(W,[1 2 4 5 6 7 8],DIM) +  TR(W,[1 2 5 6 7 8 10],DIM); %L_BE(W)=W
    TR(W,[1 2 3 4 7 8],DIM) + TR(W,[1 2 3 4 6 7 8 10],DIM) == TR(W,[1 2 3 4 6 7 8],DIM) +  TR(W,[1 2 3 4 7 8 10],DIM); %L_CE(W)=W
    TR(W,[1 2 3 4 5 6],DIM) + TR(W,[1 2 3 4 5 6 8 10],DIM) == TR(W,[1 2 3 4 5 6 8],DIM) +  TR(W,[1 2 3 4 5 6 10],DIM); %L_DE(W)=W
   

%L_ABC(W)=W 
   TR(W,[7 8 9 10],DIM) + TR(W,[2 4 7 8 9 10],DIM) +  TR(W,[4 6 7 8 9 10],DIM) +  TR(W,[2 6 7 8 9 10],DIM) ==  TR(W,[2 7 8 9 10],DIM) +  TR(W,[4 7 8 9 10],DIM) +  TR(W,[6 7 8 9 10],DIM) +  TR(W,[2 4 6 7 8 9 10],DIM);%L_{ABC}
   TR(W,[3 4 9 10],DIM) + TR(W,[2 3 4 6 9 10],DIM) +  TR(W,[2 3 4 8 9 10],DIM) +  TR(W,[3 4 6 8 9 10],DIM) ==  TR(W,[2 3 4 9 10],DIM) +  TR(W,[3 4 8 9 10],DIM) +  TR(W,[6 3 4 9 10],DIM) +  TR(W,[2 3 4 6 8 9 10],DIM);%L_{ACD}
   TR(W,[5 6 9 10],DIM) + TR(W,[2 4 5 6 9 10],DIM) +  TR(W,[4 5 6 8 9 10],DIM) +  TR(W,[2 5 6 8 9 10],DIM) ==  TR(W,[2 5 6 9 10],DIM) +  TR(W,[4 5 6 9 10],DIM) +  TR(W,[5 6 8 9 10],DIM) +  TR(W,[2 4 5 6 8 9 10],DIM);%L_{ABD}
   TR(W,[5 6 7 8],DIM) + TR(W,[2 5 6 7 8 10],DIM) +  TR(W,[4 5 6 7 8 10],DIM) +  TR(W,[2 4 5 6 7 8],DIM) ==  TR(W,[2 5 6 7 8],DIM) +  TR(W,[4 5 6 7 8],DIM) +  TR(W,[5 6 7 8 10],DIM) +  TR(W,[2 4 5 6 7 8 10],DIM);%L_{ABE}
   TR(W,[3 4 7 8],DIM) + TR(W,[2 3 4 7 8 10],DIM) +  TR(W,[3 4 6 7 8 10],DIM) +  TR(W,[2 3 4 6 7 8],DIM) ==  TR(W,[2 3 4 7 8],DIM) +  TR(W,[3 4 6 7 8],DIM) +  TR(W,[3 4 7 8 10],DIM) +  TR(W,[2 3 4 6 7 8 10],DIM);%L_{ACE}
   TR(W,[3 4 5 6],DIM) + TR(W,[2 3 4 5 6 10],DIM) +  TR(W,[3 4 5 6 8 10],DIM) +  TR(W,[2 3 4 5 6 8],DIM) ==  TR(W,[2 3 4 5 6],DIM) +  TR(W,[3 4 5 6 8],DIM) +  TR(W,[3 4 5 6 10],DIM) +  TR(W,[2 3 4 5 6 8 10],DIM);%L_{ADE}
   TR(W,[1 2 9 10],DIM) + TR(W,[1 2 4 8 9 10],DIM) +  TR(W,[1 2 6 8 9 10],DIM) +  TR(W,[1 2 4 6 9 10],DIM) ==  TR(W,[1 2 4 9 10],DIM) +  TR(W,[1 2 6 9 10],DIM) +  TR(W,[1 2 8 9 10],DIM) +  TR(W,[1 2 4 6 8 9 10],DIM);%L_{BCD}
   TR(W,[1 2 7 8],DIM) + TR(W,[1 2 4 7 8 10],DIM) +  TR(W,[1 2 6 7 8 10],DIM) +  TR(W,[1 2 4 6 7 8],DIM) ==  TR(W,[1 2 7 8 10],DIM) +  TR(W,[1 2 4 7 8],DIM) +  TR(W,[1 2 6 7 8],DIM) +  TR(W,[1 2 4 6 7 8 10],DIM);%L_{BCE}
   TR(W,[1 2 5 6],DIM) + TR(W,[1 2 4 5 6 8],DIM) +  TR(W,[1 2 4 5 6 10],DIM) +  TR(W,[1 2 5 6 8 10],DIM) ==  TR(W,[1 2 4 5 6],DIM) +  TR(W,[1 2 5 6 8],DIM) +  TR(W,[1 2 5 6 10],DIM) +  TR(W,[1 2 4 5 6 8 10],DIM);%L_{BDE}
   TR(W,[1 2 3 4],DIM) + TR(W,[1 2 3 4 8 10],DIM) +  TR(W,[1 2 3 4 6 10],DIM) +  TR(W,[1 2 3 4 6 8],DIM) ==  TR(W,[1 2 3 4 6],DIM) +  TR(W,[1 2 3 4 8],DIM) +  TR(W,[1 2 3 4 10],DIM) +  TR(W,[1 2 3 4 6 8 10],DIM);%L_{CDE}

%L_ABCD(W)=W
   TR(W,[9 10],DIM) + TR(W,[2 4 9 10],DIM) + TR(W,[4 6 9 10],DIM) + TR(W,[2 6 9 10],DIM) + TR(W,[2 8 9 10],DIM) + TR(W,[4 8 9 10],DIM) + TR(W,[6 8 9 10],DIM) + TR(W,[2 4 6 8 9 10],DIM) == TR(W,[2 9 10],DIM) + TR(W,[4 9 10],DIM) + TR(W,[6 9 10],DIM) + TR(W,[8 9 10],DIM) + TR(W,[2 4 6 9 10],DIM) + TR(W,[2 4 8 9 10],DIM) + TR(W,[2 6 8 9 10],DIM) + TR(W,[4 6 8 10],DIM); %L_{ABCD}
   TR(W,[7 8],DIM) + TR(W,[6 7 8 10],DIM) + TR(W,[4 7 8 10],DIM) + TR(W,[2 7 8 10],DIM) + TR(W,[2 4 7 8],DIM) + TR(W,[2 6 7 8],DIM) + TR(W,[4 6 7 8],DIM) + TR(W,[2 4 6 7 8 10],DIM) == TR(W,[2 7 8],DIM) + TR(W,[4 7 8],DIM) + TR(W,[6 7 8],DIM) + TR(W,[7 8 10],DIM) + TR(W,[2 4 6 7 8],DIM) + TR(W,[2 4 7 8 10],DIM) + TR(W,[2 6 7 8 10],DIM) + TR(W,[4 6 7 8 10],DIM); %L_{ABCE}
   TR(W,[5 6],DIM) + TR(W,[2 4 5 6],DIM) + TR(W,[4 5 6 10],DIM) + TR(W,[2 5 6 10],DIM) + TR(W,[2 5 6 8],DIM) + TR(W,[4 5 6 8],DIM) + TR(W,[5 6 8 10],DIM) + TR(W,[2 4 5 6 8 10],DIM) == TR(W,[2 5 6],DIM) + TR(W,[4 5 6],DIM) + TR(W,[5 6 10],DIM) + TR(W,[5 6 8],DIM) + TR(W,[2 4 5 6 8],DIM) + TR(W,[2 4 5 6 10],DIM) + TR(W,[2 5 6 8 10],DIM) + TR(W,[4 5 6 8 10],DIM); %L_{ABDE}
   TR(W,[3 4],DIM) + TR(W,[2 3 4 6],DIM) + TR(W,[3 4 8 10],DIM) + TR(W,[2 3 4 10],DIM) + TR(W,[2 3 4 8],DIM) + TR(W,[3 4 6 10],DIM) + TR(W,[3 4 6 8],DIM) + TR(W,[2 3 4 6 8 10],DIM) == TR(W,[2 3 4],DIM) + TR(W,[3 4 6],DIM) + TR(W,[3 4 10],DIM) + TR(W,[3 4 8],DIM) + TR(W,[2 3 4 6 8],DIM) + TR(W,[2 3 4 6 10],DIM) + TR(W,[2 3 4 8 10],DIM) + TR(W,[3 4 6 8 10],DIM); %L_{ACDE}
   TR(W,[1 2],DIM) + TR(W,[1 2 4 6],DIM) + TR(W,[1 2 4 8],DIM) + TR(W,[1 2 4 10],DIM) + TR(W,[1 2 6 8],DIM) + TR(W,[1 2 6 10],DIM) + TR(W,[1 2 8 10],DIM) + TR(W,[1 2 4 6 8 10],DIM) == TR(W,[1 2 10],DIM) + TR(W,[1 2 8],DIM) + TR(W,[1 2 6],DIM) + TR(W,[1 2 4],DIM) + TR(W,[1 2 4 6 8],DIM) + TR(W,[1 2 4 8 10],DIM) + TR(W,[1 2 6 8 10],DIM) + TR(W,[1 2 4 6 10],DIM); %L_{BCDE}




     W ==  -TR(W,[1 2 3 4 5 6 7 8],DIM) + TR(W,[1 2 3 4 5 6 7 8 10],DIM) -TR(W,[1 2 3 4 5 6 9 10],DIM) + TR(W,[1 2 3 4 5 6 8 9 10],DIM)-TR(W,[1 2 3 4 7 8 9 10],DIM) + TR(W,[1 2 3 4 6 7 8 9 10],DIM) -TR(W,[1 2 5 6 7 8 9 10],DIM) + TR(W,[1 2 4 5 6 7 8 9 10],DIM)  -TR(W,[3 4 5 6 7 8 9 10],DIM) + TR(W,[2 3 4 5 6 7 8 9 10],DIM)-TR(W,[5 6 7 8 9 10],DIM) - TR(W,[2 4 5 6 7 8 9 10],DIM) + TR(W,[2 5 6 7 8 9 10],DIM) + TR(W,[4 5 6 7 8 9 10],DIM) -TR(W,[1 2 3 4 9 10],DIM) - TR(W,[1 2 3 4 6 8 9 10],DIM) + TR(W,[1 2 3 4 6 9 10],DIM) + TR(W,[1 2 3 4 8 9 10],DIM) -TR(W,[1 2 5 6 9 10],DIM) - TR(W,[1 2 4 5 6 8 9 10],DIM) + TR(W,[1 2 4 5 6 9 10],DIM) +  TR(W,[1 2 5 6 8 9 10],DIM) -TR(W,[1 2 7 8 9 10],DIM) - TR(W,[1 2 4 6 7 8 9 10],DIM) + TR(W,[1 2 4 7 8 9 10],DIM) +  TR(W,[1 2 6 7 8 9 10],DIM)-TR(W,[3 4 5 6 9 10],DIM) - TR(W,[2 3 4 5 6 8 9 10],DIM) + TR(W,[2 3 4 5 6 9 10],DIM) +  TR(W,[3 4 5 6 8 9 10],DIM)-TR(W,[3 4 7 8 9 10],DIM) - TR(W,[2 3 4 6 7 8 9 10],DIM) + TR(W,[2 3 4 7 8 9 10],DIM) +  TR(W,[3 4 6 7 8 9 10],DIM)-TR(W,[3 4 5 6 7 8],DIM) - TR(W,[2 3 4 5 6 7 8 10],DIM) + TR(W,[2 3 4 5 6 7 8],DIM) +  TR(W,[3 4 5 6 7 8 10],DIM)-TR(W,[1 2 5 6 7 8],DIM) - TR(W,[1 2 4 5 6 7 8 10],DIM) + TR(W,[1 2 4 5 6 7 8],DIM) +  TR(W,[1 2 5 6 7 8 10],DIM)-TR(W,[1 2 3 4 7 8],DIM) - TR(W,[1 2 3 4 6 7 8 10],DIM) + TR(W,[1 2 3 4 6 7 8],DIM) +  TR(W,[1 2 3 4 7 8 10],DIM)-TR(W,[1 2 3 4 5 6],DIM) - TR(W,[1 2 3 4 5 6 8 10],DIM) + TR(W,[1 2 3 4 5 6 8],DIM) +  TR(W,[1 2 3 4 5 6 10],DIM)-TR(W,[7 8 9 10],DIM) - TR(W,[2 4 7 8 9 10],DIM) -  TR(W,[4 6 7 8 9 10],DIM) -  TR(W,[2 6 7 8 9 10],DIM) +  TR(W,[2 7 8 9 10],DIM) +  TR(W,[4 7 8 9 10],DIM) +  TR(W,[6 7 8 9 10],DIM) +  TR(W,[2 4 6 7 8 9 10],DIM) -TR(W,[3 4 9 10],DIM) - TR(W,[2 3 4 6 9 10],DIM) -  TR(W,[2 3 4 8 9 10],DIM) -  TR(W,[3 4 6 8 9 10],DIM) +  TR(W,[2 3 4 9 10],DIM) +  TR(W,[3 4 8 9 10],DIM) +  TR(W,[6 3 4 9 10],DIM) +  TR(W,[2 3 4 6 8 9 10],DIM) -TR(W,[5 6 9 10],DIM) - TR(W,[2 4 5 6 9 10],DIM) -  TR(W,[4 5 6 8 9 10],DIM) -  TR(W,[2 5 6 8 9 10],DIM) +  TR(W,[2 5 6 9 10],DIM) +  TR(W,[4 5 6 9 10],DIM) +  TR(W,[5 6 8 9 10],DIM) +  TR(W,[2 4 5 6 8 9 10],DIM) -TR(W,[5 6 7 8],DIM) - TR(W,[2 5 6 7 8 10],DIM) -  TR(W,[4 5 6 7 8 10],DIM) -  TR(W,[2 4 5 6 7 8],DIM) +  TR(W,[2 5 6 7 8],DIM) +  TR(W,[4 5 6 7 8],DIM) +  TR(W,[5 6 7 8 10],DIM) +  TR(W,[2 4 5 6 7 8 10],DIM) -TR(W,[3 4 7 8],DIM) - TR(W,[2 3 4 7 8 10],DIM) -  TR(W,[3 4 6 7 8 10],DIM) -  TR(W,[2 3 4 6 7 8],DIM) +  TR(W,[2 3 4 7 8],DIM) +  TR(W,[3 4 6 7 8],DIM) +  TR(W,[3 4 7 8 10],DIM) +  TR(W,[2 3 4 6 7 8 10],DIM) -TR(W,[3 4 5 6],DIM) - TR(W,[2 3 4 5 6 10],DIM) -  TR(W,[3 4 5 6 8 10],DIM) -  TR(W,[2 3 4 5 6 8],DIM) +  TR(W,[2 3 4 5 6],DIM) +  TR(W,[3 4 5 6 8],DIM) +  TR(W,[3 4 5 6 10],DIM) +  TR(W,[2 3 4 5 6 8 10],DIM)-TR(W,[1 2 9 10],DIM) - TR(W,[1 2 4 8 9 10],DIM) -  TR(W,[1 2 6 8 9 10],DIM) -  TR(W,[1 2 4 6 9 10],DIM) +  TR(W,[1 2 4 9 10],DIM) +  TR(W,[1 2 6 9 10],DIM) +  TR(W,[1 2 8 9 10],DIM) +  TR(W,[1 2 4 6 8 9 10],DIM)-TR(W,[1 2 7 8],DIM) - TR(W,[1 2 4 7 8 10],DIM) -  TR(W,[1 2 6 7 8 10],DIM) -  TR(W,[1 2 4 6 7 8],DIM) +  TR(W,[1 2 7 8 10],DIM) +  TR(W,[1 2 4 7 8],DIM) +  TR(W,[1 2 6 7 8],DIM) +  TR(W,[1 2 4 6 7 8 10],DIM)-TR(W,[1 2 5 6],DIM) - TR(W,[1 2 4 5 6 8],DIM) -  TR(W,[1 2 4 5 6 10],DIM) -  TR(W,[1 2 5 6 8 10],DIM) +  TR(W,[1 2 4 5 6],DIM) +  TR(W,[1 2 5 6 8],DIM) +  TR(W,[1 2 5 6 10],DIM) +  TR(W,[1 2 4 5 6 8 10],DIM)-TR(W,[1 2 3 4],DIM) - TR(W,[1 2 3 4 8 10],DIM) -  TR(W,[1 2 3 4 6 10],DIM) -  TR(W,[1 2 3 4 6 8],DIM) +  TR(W,[1 2 3 4 6],DIM) +  TR(W,[1 2 3 4 8],DIM) +  TR(W,[1 2 3 4 10],DIM) +  TR(W,[1 2 3 4 6 8 10],DIM)-TR(W,[9 10],DIM) - TR(W,[2 4 9 10],DIM) - TR(W,[4 6 9 10],DIM) - TR(W,[2 6 9 10],DIM) - TR(W,[2 8 9 10],DIM) - TR(W,[4 8 9 10],DIM) - TR(W,[6 8 9 10],DIM) - TR(W,[2 4 6 8 9 10],DIM) + TR(W,[2 9 10],DIM) + TR(W,[4 9 10],DIM) + TR(W,[6 9 10],DIM) + TR(W,[8 9 10],DIM) + TR(W,[2 4 6 9 10],DIM) + TR(W,[2 4 8 9 10],DIM) + TR(W,[2 6 8 9 10],DIM) + TR(W,[4 6 8 10],DIM) -TR(W,[7 8],DIM) - TR(W,[6 7 8 10],DIM) - TR(W,[4 7 8 10],DIM) - TR(W,[2 7 8 10],DIM) - TR(W,[2 4 7 8],DIM) - TR(W,[2 6 7 8],DIM) - TR(W,[4 6 7 8],DIM) - TR(W,[2 4 6 7 8 10],DIM) + TR(W,[2 7 8],DIM) + TR(W,[4 7 8],DIM) + TR(W,[6 7 8],DIM) + TR(W,[7 8 10],DIM) + TR(W,[2 4 6 7 8],DIM) + TR(W,[2 4 7 8 10],DIM) + TR(W,[2 6 7 8 10],DIM) + TR(W,[4 6 7 8 10],DIM) -TR(W,[5 6],DIM) - TR(W,[2 4 5 6],DIM) - TR(W,[4 5 6 10],DIM) - TR(W,[2 5 6 10],DIM) - TR(W,[2 5 6 8],DIM) - TR(W,[4 5 6 8],DIM) - TR(W,[5 6 8 10],DIM) - TR(W,[2 4 5 6 8 10],DIM) + TR(W,[2 5 6],DIM) + TR(W,[4 5 6],DIM) + TR(W,[5 6 10],DIM) + TR(W,[5 6 8],DIM) + TR(W,[2 4 5 6 8],DIM) + TR(W,[2 4 5 6 10],DIM) + TR(W,[2 5 6 8 10],DIM) + TR(W,[4 5 6 8 10],DIM) -TR(W,[3 4],DIM) - TR(W,[2 3 4 6],DIM) - TR(W,[3 4 8 10],DIM) - TR(W,[2 3 4 10],DIM) - TR(W,[2 3 4 8],DIM) - TR(W,[3 4 6 10],DIM) - TR(W,[3 4 6 8],DIM) - TR(W,[2 3 4 6 8 10],DIM) + TR(W,[2 3 4],DIM) + TR(W,[3 4 6],DIM) + TR(W,[3 4 10],DIM) + TR(W,[3 4 8],DIM) + TR(W,[2 3 4 6 8],DIM) + TR(W,[2 3 4 6 10],DIM) + TR(W,[2 3 4 8 10],DIM) + TR(W,[3 4 6 8 10],DIM)-TR(W,[1 2],DIM) + TR(W,[1 2 4 6],DIM) + TR(W,[1 2 4 8],DIM) + TR(W,[1 2 4 10],DIM) + TR(W,[1 2 6 8],DIM) + TR(W,[1 2 6 10],DIM) + TR(W,[1 2 8 10],DIM) + TR(W,[1 2 4 6 8 10],DIM) + TR(W,[1 2 10],DIM) + TR(W,[1 2 8],DIM) + TR(W,[1 2 6],DIM) + TR(W,[1 2 4],DIM) + TR(W,[1 2 4 6 8],DIM) + TR(W,[1 2 4 8 10],DIM) + TR(W,[1 2 6 8 10],DIM) + TR(W,[1 2 4 6 10],DIM);

       
     trace(W) == dOut^k;

     
%%% Conditions for k=2 and k=3
     %  TR(W,[1 2],DIM) == TR(W,[1 2 4],DIM);
     % TR(W,[3 4],DIM) == TR(W,[2 3 4],DIM);
     % W == TR(W,2,DIM) + TR(W,4,DIM) - TR(W,[2 4],DIM);
     % trace(W) == dOut^2;

     % TR(W,[1 2 3 4],DIM) == TR(W,[1 2 3 4 6],DIM);
     % TR(W,[3 4 5 6],DIM) == TR(W,[2 3 4 5 6],DIM);
     % TR(W,[1 2 5 6],DIM) == TR(W,[1 2 4 5 6],DIM);
     % TR(W,[1 2],DIM) + TR(W,[1 2 4 6],DIM) ==  TR(W,[1 2 4],DIM) +  TR(W,[1 2 6],DIM);
     % TR(W,[3 4],DIM) + TR(W,[2 3 4 6],DIM) ==  TR(W,[2 3 4],DIM) +  TR(W,[3 4 6],DIM);
     % TR(W,[5 6],DIM) + TR(W,[2 4 5 6],DIM) ==  TR(W,[4 5 6],DIM) +  TR(W,[2 5 6],DIM);
     % W == TR(W,[2 4 6],DIM) + TR(W,2,DIM) + TR(W,4,DIM) + TR(W,6,DIM) - TR(W,[2 4],DIM) - TR(W,[2 6],DIM) - TR(W,[4 6],DIM);
     % trace(W) == dOut^3;
  
     maximise pS;  
cvx_end   
%%%%%%%%%%%%%%%%%%%%% ERROR MESSAGE %%%%%%%%%%%%%%%%%%%%%
    otherwise
        'ERROR!!'
        'Set protocol equals 1 for PAR, 2 for SEQ, and 4 for GEN'
        pause
end

end