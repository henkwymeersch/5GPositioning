function [indices,h]=DCSSOMP(Y,A,L)
% This function implements the paper
% Sarvotham, S., Baron, D., Wakin, M., Duarte, M. F., & Baraniuk, R. G. (2005, November). 
% "Distributed compressed sensing of jointly sparse signals." 
% In Asilomar conference on signals, systems, and computers (pp. 1537-1541)
% (numbering refers to numbers in the paper, section IV.B)
% 
% usage: [indices,h]=DCSSOMP(Y,A,L)
% Y = input (one column vector for each channel / subcarrier)
% A = sensing matrix for each channel
% L = the sparsity level
% 
% output: 
% indices = index of L chosen columns of A
% h = the recovered coefficients
%
% code version 1.0, September 2018, Henk Wymeersch & Gonzalo Seco Granados, henkw@chalmers.se

    K=size(Y,2);        % number of channels
    N=size(Y,1);        % observations per channel
    M=size(A,2);        % size of sparse vector (M>>N)
    if (L<=0)
        L=N;
    end    
    % 1. initialization 
    R=Y;                % residual
    psi=zeros(N,L,K);
    indices=zeros(1,L);
    columns=zeros(N,L,K);
    betamatrix=zeros(L,K);  
    
    for counter=1:L  
        % 2. find maximum correlation between residual and columns of A
        cost=zeros(1,M);
        for m=1:M        
            for k=1:K
                cost(m)=cost(m)+abs(A(:,m,k)'*R(:,k))/norm(A(:,m,k));
            end
        end
        [maxval,maxi]=max(cost);   
        indices(counter)=maxi;   
        
        for k=1:K
            % 3. orthogonalize
            columns(:,counter,k)=A(:,maxi,k);
            omega=A(:,maxi,k);                        
            psi(:,counter,k)=omega;
            for counter2=1:counter-1                
                psi(:,counter,k)=psi(:,counter,k)-(psi(:,counter2,k)'*omega)*psi(:,counter2,k)/(norm(psi(:,counter2,k)))^2; 
            end  
            
            % 4. update coefficients and residual
            beta=psi(:,counter,k)'*R(:,k)/(norm(psi(:,counter,k)))^2;
            betamatrix(counter,k)=beta;
            R(:,k)=R(:,k)-beta*psi(:,counter,k);
        end                  
    end
    
    % 6. Deorthogonalize 
    h=zeros(L,K);
    for k=1:K
        [Q,Rqr]=qr(columns(:,:,k),0);    
        h(:,k)=inv(Rqr)*Q'*psi(:,:,k)*betamatrix(:,k);    
    end       
