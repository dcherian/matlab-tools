function [dPdx,dPdy] = roms_prsgrd32(rho,zrmat,zwmat,Hz)

%svn $Id: prsgrd.F 645 2013-01-22 23:21:54Z arango $
%================================================== Hernan G. Arango ===
%  Copyright (c) 2002-2013 The ROMS/TOMS Group                         %
%    Licensed under a MIT/X style license                              %
%    See License_ROMS.txt                                              %
%=======================================================================
%                                                                      %
%  This routine computes the baroclinic hydrostatic pressure gradient  %
%  term.                                                               %
%                                                                      %
%=======================================================================
%
%
%svn $Id: prsgrd32.h 645 2013-01-22 23:21:54Z arango $
%***********************************************************************
%  Copyright (c) 2002-2013 The ROMS/TOMS Group                         %
%    Licensed under a MIT/X style license                              %
%    See License_ROMS.txt                           Hernan G. Arango   %
%****************************************** Alexander F. Shchepetkin ***
%                                                                      %
%  This subroutine evaluates the nonlinear  baroclinic,  hydrostatic   %
%  pressure gradient term using a  nonconservative  Density-Jacobian   %
%  scheme,  based on  cubic polynomial fits for  "rho" and  "z_r" as   %
%  functions of nondimensional coordinates (XI,ETA,s), that is,  its   %
%  respective array indices. The  cubic polynomials  are monotonized   %
%  by using  harmonic mean instead of linear averages to interpolate   %
%  slopes. This scheme retains exact anti-symmetry:                    %
%                                                                      %
%        J(rho,z_r)=-J(z_r,rho).                                       %
%                                                                      %
%  If parameter OneFifth (below) is set to zero,  the scheme becomes   %
%  identical to standard Jacobian.                                     %
%                                                                      %
%  Reference:                                                          %
%                                                                      %
%    Shchepetkin A.F and J.C. McWilliams, 2003:  A method for          %
%      computing horizontal pressure gradient force in an ocean        %
%      model with non-aligned vertical coordinate, JGR, 108,           %
%      1-34.                                                           %
%                                                                      %
%***********************************************************************
%
%  Set horizontal starting and ending indices for automatic private
%  storage arrays.
%
%       IminS=BOUNDS(ng)%Istr(tile)-3
%       ImaxS=BOUNDS(ng)%Iend(tile)+3
%       JminS=BOUNDS(ng)%Jstr(tile)-3
%       JmaxS=BOUNDS(ng)%Jend(tile)+3
% %
% %  Determine array lower and upper bounds in the I- and J-directions.
% %
%       LBi=BOUNDS(ng)%LBi(tile)
%       UBi=BOUNDS(ng)%UBi(tile)
%       LBj=BOUNDS(ng)%LBj(tile)
%       UBj=BOUNDS(ng)%UBj(tile)
% %
% %  Set array lower and upper bounds for MIN(I,J) directions and
% %  MAX(I,J) directions.
% %
%       LBij=BOUNDS(ng)%LBij
%       UBij=BOUNDS(ng)%UBij

%%  Imported variable declarations.
%
%       integer, intent(in) :: ng, tile
%       integer, intent(in) :: LBi, UBi, LBj, UBj
%       integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
%       integer, intent(in) :: nrhs
%       real(r8), intent(in) :: om_v(LBi:,LBj:)
%       real(r8), intent(in) :: on_u(LBi:,LBj:)
%       real(r8), intent(in) :: Hz(LBi:,LBj:,:)
%       real(r8), intent(in) :: z_r(LBi:,LBj:,:)
%       real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
%       real(r8), intent(in) :: rho(LBi:,LBj:,:)
%       real(r8), intent(inout) :: ru(LBi:,LBj:,0:,:)
%       real(r8), intent(inout) :: rv(LBi:,LBj:,0:,:)
%
%  Local variable declarations.
%
      
%       real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: P
%       real(r8), dimension(IminS:ImaxS,0:N(ng)) :: dR
%       real(r8), dimension(IminS:ImaxS,0:N(ng)) :: dZ
%       real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FC
%       real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: aux
%       real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: dRx
%       real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: dZx
% %
% %
% %-----------------------------------------------------------------------
% %  Set lower and upper tile bounds and staggered variables bounds for
% %  this horizontal domain partition.  Notice that if tile=-1, it will
% %  set the values for the global grid.
% %-----------------------------------------------------------------------
% %
%       integer :: Istr, IstrB, IstrP, IstrR, IstrT, IstrM, IstrU
%       integer :: Iend, IendB, IendP, IendR, IendT
%       integer :: Jstr, JstrB, JstrP, JstrR, JstrT, JstrM, JstrV
%       integer :: Jend, JendB, JendP, JendR, JendT
%       integer :: Istrm3, Istrm2, Istrm1, IstrUm2, IstrUm1
%       integer :: Iendp1, Iendp2, Iendp2i, Iendp3
%       integer :: Jstrm3, Jstrm2, Jstrm1, JstrVm2, JstrVm1
%       integer :: Jendp1, Jendp2, Jendp2i, Jendp3
% %
%       Istr   =BOUNDS(ng) % Istr   (tile)
%       IstrB  =BOUNDS(ng) % IstrB  (tile)
%       IstrM  =BOUNDS(ng) % IstrM  (tile)
%       IstrP  =BOUNDS(ng) % IstrP  (tile)
%       IstrR  =BOUNDS(ng) % IstrR  (tile)
%       IstrT  =BOUNDS(ng) % IstrT  (tile)
%       IstrU  =BOUNDS(ng) % IstrU  (tile)
%       Iend   =BOUNDS(ng) % Iend   (tile)
%       IendB  =BOUNDS(ng) % IendB  (tile)
%       IendP  =BOUNDS(ng) % IendP  (tile)
%       IendR  =BOUNDS(ng) % IendR  (tile)
%       IendT  =BOUNDS(ng) % IendT  (tile)
%       Jstr   =BOUNDS(ng) % Jstr   (tile)
%       JstrB  =BOUNDS(ng) % JstrB  (tile)
%       JstrM  =BOUNDS(ng) % JstrM  (tile)
%       JstrP  =BOUNDS(ng) % JstrP  (tile)
%       JstrR  =BOUNDS(ng) % JstrR  (tile)
%       JstrT  =BOUNDS(ng) % JstrT  (tile)
%       JstrV  =BOUNDS(ng) % JstrV  (tile)
%       Jend   =BOUNDS(ng) % Jend   (tile)
%       JendB  =BOUNDS(ng) % JendB  (tile)
%       JendP  =BOUNDS(ng) % JendP  (tile)
%       JendR  =BOUNDS(ng) % JendR  (tile)
%       JendT  =BOUNDS(ng) % JendT  (tile)
% %
%       Istrm3 =BOUNDS(ng) % Istrm3 (tile)            % Istr-3
%       Istrm2 =BOUNDS(ng) % Istrm2 (tile)            % Istr-2
%       Istrm1 =BOUNDS(ng) % Istrm1 (tile)            % Istr-1
%       IstrUm2=BOUNDS(ng) % IstrUm2(tile)            % IstrU-2
%       IstrUm1=BOUNDS(ng) % IstrUm1(tile)            % IstrU-1
%       Iendp1 =BOUNDS(ng) % Iendp1 (tile)            % Iend+1
%       Iendp2 =BOUNDS(ng) % Iendp2 (tile)            % Iend+2
%       Iendp2i=BOUNDS(ng) % Iendp2i(tile)            % Iend+2 interior
%       Iendp3 =BOUNDS(ng) % Iendp3 (tile)            % Iend+3
%       Jstrm3 =BOUNDS(ng) % Jstrm3 (tile)            % Jstr-3
%       Jstrm2 =BOUNDS(ng) % Jstrm2 (tile)            % Jstr-2
%       Jstrm1 =BOUNDS(ng) % Jstrm1 (tile)            % Jstr-1
%       JstrVm2=BOUNDS(ng) % JstrVm2(tile)            % JstrV-2
%       JstrVm1=BOUNDS(ng) % JstrVm1(tile)            % JstrV-1
%       Jendp1 =BOUNDS(ng) % Jendp1 (tile)            % Jend+1
%       Jendp2 =BOUNDS(ng) % Jendp2 (tile)            % Jend+2
%       Jendp2i=BOUNDS(ng) % Jendp2i(tile)            % Jend+2 interior
%       Jendp3 =BOUNDS(ng) % Jendp3 (tile)            % Jend+3
%

%% preliminary code

%-----------------------------------------------------------------------
%  Preliminary step (same for XI- and ETA-components:
%-----------------------------------------------------------------------
% Required constants
    g = 9.81;
    rho0 = 1025;
    
    GRho = g/rho0;
    GRho0 = 1000.0*GRho;
    HalfGRho=0.5 *GRho;

    OneFifth = 0.2;
    OneTwelfth = 1.0/12.0;
    eps = 1.0E-10;
    
% allocate variables    
    [L,M,N] = size(zrmat);
    
    dR = nan(size(zwmat));
    dZ = nan(size(zwmat));
    P = nan(size(rho));

% All steps follow Shchepetkin & McWilliams (2003)
% Step 1: First set dR and dZ (provisional)
    dR(:,:,2:end-1) = diff(rho,1,3);
    dZ(:,:,2:end-1) = diff(zrmat,1,3);

    dR(:,:,end) = dR(:,:,end-1);
    dZ(:,:,end) = dZ(:,:,end-1);

    dR(:,:,1) = dR(:,:,2);
    dZ(:,:,1) = dZ(:,:,2);
      
% Step 2: Compute harmonic averages
    for k=2:N
        dR1(:,:,k) = dR(:,:,k) .* dR(:,:,k-1);
        dZ1(:,:,k) = dZ(:,:,k) .* dZ(:,:,k-1);
    end
    dZ1 = dZ1 ./ avg1(dZ,3);
    dR1(dR1 < eps) = 0;
    dR1 = dR1 ./ avg1(dR,3);
    
    clear dR dZ
    dR = dR1; dZ = dZ1;
    clear dR1 dZ1
    
% Step 3: Apply boundary conditions at the top

% Step 4: top most grid box
    cff1=1 ./ (zrmat(:,:,end)-zrmat(:,:,end-1));
    cff2=0.5*(rho(:,:,end)-rho(:,:,end-1)) .*   ...              
            (zwmat(:,:,end)-zrmat(:,:,end)) .*cff1;
        
    P(:,:,end) = GRho0*zwmat(:,:,end)+           ...                  
                      GRho*(rho(:,:,end)+cff2) .*  ...                    
                      (zwmat(:,:,end)-zrmat(:,:,end));
    
    for k= N-1:-1:1
        P(:,:,k)=P(:,:,k+1)+ ...                                       
                    HalfGRho*((rho(:,:,k+1)+rho(:,:,k)) .* ...              
                              (zrmat(:,:,k+1)-zrmat(:,:,k))-    ...           
                              OneFifth*                        ...        
                              ((dR(:,:,k+1)-dR(:,:,k)) .*               ...     
                               (zrmat(:,:,k+1)-zrmat(:,:,k)- ...              
                                OneTwelfth*                 ...           
                                (dZ(:,:,k+1)+dZ(:,:,k)))-          ...        
                               (dZ(:,:,k+1)-dZ(:,:,k)) .*               ...     
                               (rho(:,:,k+1)-rho(:,:,k)-             ...  
                                OneTwelfth*              ...              
                                (dR(:,:,k+1)+dR(:,:,k)))));
    end
    
%%      for j=JstrV-1:Jend
%         for k=N(ng):-1:1
%           for i=IstrU-1:Iend
%             cff=2.0*dR(i,k)*dR(i,k-1);
%             if (cff > eps) 
%               dR(i,k)=cff/(dR(i,k)+dR(i,k-1));
%             else
%               dR(i,k)=0.0;
%             end
%             dZ(i,k)=2.0*dZ(i,k)*dZ(i,k-1)/(dZ(i,k)+dZ(i,k-1));
%           end
%         end
        
%         for i=IstrU-1:Iend
%           cff1=1.0/(z_r(i,j,N(ng))-z_r(i,j,N(ng)-1));
%           cff2=0.5*(rho(i,j,N(ng))-rho(i,j,N(ng)-1))*   ...              
%               (z_w(i,j,N(ng))-z_r(i,j,N(ng)))*cff1;
%           P(i,j,N(ng))=GRho0*z_w(i,j,N(ng))+          ...                  
%                       GRho*(rho(i,j,N(ng))+cff2)*  ...                    
%                       (z_w(i,j,N(ng))-z_r(i,j,N(ng)));
%         end
%         for k=N(ng)-1:-1:1
%           for i=IstrU-1:Iend
%             P(i,j,k)=P(i,j,k+1)+ ...                                       
%                     HalfGRho*((rho(i,j,k+1)+rho(i,j,k))* ...              
%                               (z_r(i,j,k+1)-z_r(i,j,k))-    ...           
%                               OneFifth*                        ...        
%                               ((dR(i,k+1)-dR(i,k))*               ...     
%                                (z_r(i,j,k+1)-z_r(i,j,k)- ...              
%                                 OneTwelfth*                 ...           
%                                 (dZ(i,k+1)+dZ(i,k)))-          ...        
%                                (dZ(i,k+1)-dZ(i,k))*               ...     
%                                (rho(i,j,k+1)-rho(i,j,k)-             ...  
%                                 OneTwelfth*              ...              
%                                 (dR(i,k+1)+dR(i,k)))));
%           end
%         end
%       end
%
%-----------------------------------------------------------------------
%%  Compute XI-component pressure gradient term.
%-----------------------------------------------------------------------
%

dZx = diff(zrmat,1,1);
dRx = diff(rho,1,1);

for i=2:L-1
    dRx1(i,:,:) = dRx(i,:,:) .* dRx(i+1,:,:);
    dZx1(i,:,:) = dZx(i,:,:) .* dZx(i+1,:,:);
end

dZx1 = dZx1 ./ avg1(dZx,1);
dRx1 = dRx1 ./ avg1(dRx,1);

dZx(dZx < eps) = 0;
dRx(dRx < eps) = 0;

for i=1:L
    dPdx(i,:,:)=on_u(i,:,:)*0.5*        ...                    
                          (Hz(i,:,:)+Hz(i-1,:,:))*  ...                    
                          (P(i-1,:,:)-P(i,:,:)-    ...                    
                           HalfGRho*              ...                     
                           ((rho(i,:,:)+rho(i-1,:,:))* ...           
                            (z_r(i,:,:)-z_r(i-1,:,:))-    ...             
                             OneFifth*                       ...          
                             ((dRx(i,:,:)-dRx(i-1,:,:))*            ...       
                              (z_r(i,:,:)-z_r(i-1,:,:)-            ...    
                               OneTwelfth*            ...                 
                               (dZx(i,:,:)+dZx(i-1,:,:)))-   ...              
                              (dZx(i,:,:)-dZx(i-1,:,:))*        ...           
                              (rho(i,:,:)-rho(i-1,:,:)-        ...        
                               OneTwelfth*                        ...     
                               (dRx(i,:,:)+dRx(i-1,:,:))))));
end


%       for k=N:-1:1
%         for j=Jstr,Jend
%           for i=IstrU-1,Iend
%             cff=2.0*aux(i,j)*aux(i+1,j);
%             if (cff > eps) 
%               cff1=1.0/(aux(i,j)+aux(i+1,j));
%               dZx(i,j)=cff*cff1;
%             else
%               dZx(i,j)=0.0;
%             end
%             cff1=2.0*FC(i,j)*FC(i+1,j);
%             if (cff1 > eps) 
%               cff2=1.0/(FC(i,j)+FC(i+1,j));
%               dRx(i,j)=cff1*cff2;
%             else
%               dRx(i,j)=0.0;
%             end
%           end
%         end

%         for j=Jstr:Jend
%           for i=IstrU:Iend
%             ru(i,j,k,nrhs)=on_u(i,:)*0.5*        ...                    
%                           (Hz(i,:,:)+Hz(i-1,j,k))*  ...                    
%                           (P(i-1,j,k)-P(i,j,k)-    ...                    
%                            HalfGRho*              ...                     
%                            ((rho(i,j,k)+rho(i-1,j,k))* ...           
%                             (z_r(i,j,k)-z_r(i-1,j,k))-    ...             
%                              OneFifth*                       ...          
%                              ((dRx(i,j)-dRx(i-1,j))*            ...       
%                               (z_r(i,j,k)-z_r(i-1,j,k)-            ...    
%                                OneTwelfth*            ...                 
%                                (dZx(i,j)+dZx(i-1,j)))-   ...              
%                               (dZx(i,j)-dZx(i-1,j))*        ...           
%                               (rho(i,j,k)-rho(i-1,j,k)-        ...        
%                                OneTwelfth*                        ...     
%                                (dRx(i,j)+dRx(i-1,j))))));
%           end
%         end
% %       end
% %
%-----------------------------------------------------------------------
%%  ETA-component pressure gradient term.
%-----------------------------------------------------------------------
%
%       for k=N(ng),1,-1
%         for j=JstrV-1,Jend+1
%           for i=Istr,Iend
%             aux(i,j)=z_r(i,j,k)-z_r(i,j-1,k)
%             FC(i,j)=rho(i,j,k)-rho(i,j-1,k)
%           end
%         end
% %
%         for j=JstrV-1,Jend
%           for i=Istr,Iend
%             cff=2.0*aux(i,j)*aux(i,j+1);
%             if (cff > eps) 
%               cff1=1.0/(aux(i,j)+aux(i,j+1));
%               dZx(i,j)=cff*cff1;
%             else
%               dZx(i,j)=0.0;
%             end
%             cff1=2.0*FC(i,j)*FC(i,j+1);
%             if (cff1 > eps) 
%               cff2=1.0/(FC(i,j)+FC(i,j+1));
%               dRx(i,j)=cff1*cff2;
%             else
%               dRx(i,j)=0.0;
%             end
%           end
%         end
% %
%         for j=JstrV,Jend
%           for i=Istr,Iend
%             rv(i,j,k,nrhs)=om_v(i,j)*0.5*   ...                        
%                           (Hz(i,j,k)+Hz(i,j-1,k))*  ...                   
%                           (P(i,j-1,k)-P(i,j,k)-        ...                
%                            HalfGRho*                      ...             
%                            ((rho(i,j,k)+rho(i,j-1,k))*       ...          
%                             (z_r(i,j,k)-z_r(i,j-1,k))-          ...       
%                              OneFifth*                             ...    
%                              ((dRx(i,j)-dRx(i,j-1))*                  ... 
%                               (z_r(i,j,k)-z_r(i,j-1,k)-    ...            
%                                OneTwelfth*                    ...         
%                                (dZx(i,j)+dZx(i,j-1)))-           ...      
%                               (dZx(i,j)-dZx(i,j-1))*                ...   
%                               (rho(i,j,k)-rho(i,j-1,k)-                ...
%                                OneTwelfth*              ...               
%                                (dRx(i,j)+dRx(i,j-1))))));
%           end
%         end
%       end
