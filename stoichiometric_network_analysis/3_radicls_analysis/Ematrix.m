% DETERMINATION OF THE EXTREME CURRENTS IN THE SNA THEORY
% Guy Schmitz, 2008
% ---------------------------------------------------------
% The program Ematrix.m find all positive solutions Ei of the equations SC*Ei = 0 where Ei
% is an extreme current and SC the matrix of the stoichiometric coefficients.
% The matrix SC of the studied model must be provided hereafter by the user of the program.
% The given example is used in the paper "Stoichiometric network analysis and associated
% dimensionless kinetic equations" by Guy Schmitz, Ljiljana Kolar-Anic, Slobodan Anic and
% Zeljko Cupic, 2008.
%
% MAIN VARIABLES
% Nr: number of reactions.
% Ra: rank of SC. This program assumes that SC is not singular.
% A(Nr,Nr): a matrix built from SC with three parts. The first lines are the lines of SC. The
% following lines (noted A1) set Nr-Ra-1 components of Ei to zero, the last line (noted A2)
% sets one component to one to avoid the trivial solution.
% B(Nr): A vector of zeros except the last element equal to one.
%
% In order to find all the possible Ei, the program builds all the possible matrix A(Nr,Nr) and
% solves the equations A*Ei=B. The main loop computes A using all the possible combinations of Nr-Ra-1
% zeros starting from all binary numbers (Bin) having Nr digits and Nr-Ra-1 ones.
% ------------------------------------------------------------------------------------------
clear
%
% % R n°    r1 r2 r3 r4 r5 r8 r9 r16
% SC(1,:)=[-1  0  0  1  0  0  0  0 ];   % Y
% SC(2,:)=[ 1 -1  0  0  0  0  0  0 ];   % M
% SC(3,:)=[ 0  1  0 -1  0  1  0  0 ];   % D
% SC(4,:)=[ 0 -1  2 -1  0  0  1  0 ];   % Z
% SC(5,:)=[ 1  0 -1  0  0 -1 -1 -2 ];   % W
% SC(6,:)=[ 0  1 -1  0 -2  0  0  0 ];   % H
load("test.mat")
SC = [-1  0  0  1  0  0  0  0; 
       1 -1  0  0  0  0  0  0;
       0  1  0 -1  0  1  0  0;
       0 -1  2 -1  0  0  1  0;
       1  0 -1  0  0 -1 -1 -2;
       0  1 -1  0 -2  0  0  0];
SC =  double(S); 
%SC(1,:)=[-1  0  0  1  0  0  0  0  0 -1];   % Y
%SC(2,:)=[ 1 -1  0  0  0  0  0  0  1  0];   % M
%SC(3,:)=[ 0  1  0 -1  0  1  0  0  0  1];   % D
%SC(4,:)=[ 0 -1  2 -1  0  0  1  0 -1  1];   % Z
%SC(5,:)=[ 1  0 -1  0  0 -1 -1 -2  0  0];   % W
%SC(6,:)=[ 0  1 -1  0 -2  0  0  0  0  0];   % H
% SC(7,:)=[0 0  0  0 -1];
% r9 es Fe2+ + NO -> M
Tiny=1e-10; % A very small number larger than the numerical errors.
Nr=size(SC,2);
Ra=rank(SC);
if Ra<size(SC,1)
   error('Matrix SC is singular')
end
B=zeros(Nr,1); B(Nr)=1;
N1=Nr-Ra-1;
Max=2^Nr-2^(Nr-N1); % Largest used binary number.
Nsol=1;
%
% Start of the main loop.
%
for K=1:Max
   Bin=dec2bin(K,Nr);
if sum(str2num(Bin(:)))==N1
   A1=zeros(N1,Nr);
   LA=1;
   for J=1:Nr
	  if Bin(J)=='1'
		A1(LA,J)=1; LA=LA+1;
	  end
	end
% Start searching one solution for a given A1
Found=0;
A2=zeros(1,Nr);
	for ColA=1:Nr
   	A2(ColA)=1;
      A=[SC;A1;A2];
      if and(Found==0,rank(A)==Nr)
         Ej=A\B;
         Found=1;
         if min(Ej)>(-Tiny)
            Add=1;
            if Nsol>1
               for ColE=1:Nsol-1
                  if (Ej-E(:,ColE))<(Tiny*ones(Nr,1))
                     Add=0; % The solution is equal to a former one.
                  end
               end
            end
            if Add==1
            	E(:,Nsol)=Ej; % Adds a new non-negative solution to E
               % disp(['K= ' num2str(K) '  Solution ' num2str(Nsol)]);
               Nsol=Nsol+1;
            end
         end
       end
  	end
% End searching one solution
end; end
%
% End of the main loop
%
% Scaling the extreme currents
for I=1:Nsol-1
   Min=1;
   for J=1:Nr
      if and(E(J,I)>Tiny, E(J,I)<Min)
         Min=E(J,I);
      end
   end
   E(:,I)=E(:,I)./Min;
end
disp(round(E))

