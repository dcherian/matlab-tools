function [dPdx,dPdy] = roms_prsgrd32v2(rho,xrmat,yrmat,z_r,z_w)

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


%% Imported variable declarations.
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
%       integer :: i, j, k
%       real(r8), parameter :: OneFifth = 0.2
%       real(r8), parameter :: OneTwelfth = 1.0/12.0
%       real(r8), parameter :: eps = 1.0E-10
%       real(r8) :: GRho, GRho0,  HalfGRho
%       real(r8) :: cff, cff1, cff2
%       real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N) :: P
%       real(r8), dimension(IminS:ImaxS,0:N) :: dR
%       real(r8), dimension(IminS:ImaxS,0:N) :: dZ
%       real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FC
%       real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: aux
%       real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: dRx
%       real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: dZx
%
%
%-----------------------------------------------------------------------
%  Set lower and upper tile bounds and staggered variables bounds for
%  this horizontal domain partition.  Notice that if tile=-1, it will
%  set the values for the global grid.
%-----------------------------------------------------------------------
%
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

%% 
%-----------------------------------------------------------------------
%  Preliminary step (same for XI- and ETA-components:
%-----------------------------------------------------------------------
%

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
    [L,M,N] = size(z_r);
    
    % calculate required grid variables
    Hz = diff(z_w,1,3);
    on_u = diff(yrmat,1,2);
    om_v = diff(xrmat,1,1);
    
%     dR = nan(size(z_w));
%     dZ = nan(size(z_w));
%     P = nan(size(rho));
    
    % Loop indices - pg. 83 of ROMS Manual (2012)
    Istr = 2;
    IstrU = 3;
    Iend = L-1;
    
    Jstr = 2;
    JstrV = 3;
    Jend = M-1;
        
    % Common steps
      for j=JstrV-1:Jend
        for k=1:N-1
          for i=IstrU-1:Iend
            dR(i,k)=rho(i,j,k+1)-rho(i,j,k);
            dZ(i,k)=z_r(i,j,k+1)-z_r(i,j,k);
          end
        end
        for i=IstrU-1:Iend
          dR(i,N)=dR(i,N-1);
          dZ(i,N)=dZ(i,N-1);
          dR(i,1)=dR(i,2);
          dZ(i,1)=dZ(i,2);
        end
        for k=N:-1:2
          for i=IstrU-1:Iend
            cff=2.0*dR(i,k)*dR(i,k-1);
            if (cff > eps) 
              dR(i,k)=cff/(dR(i,k)+dR(i,k-1));
            else
              dR(i,k)=0.0;
            end
            dZ(i,k)=2.0*dZ(i,k)*dZ(i,k-1)/(dZ(i,k)+dZ(i,k-1));
          end
        end
        for i=IstrU-1:Iend
          cff1=1.0/(z_r(i,j,N)-z_r(i,j,N-1));
          cff2=0.5*(rho(i,j,N)-rho(i,j,N-1)) .*               ...
              (z_w(i,j,N)-z_r(i,j,N))*cff1;                   ...
         P(i,j,N)=GRho0*z_w(i,j,N)+                            ...
                      GRho*(rho(i,j,N)+cff2) .*                      ...
                      (z_w(i,j,N)-z_r(i,j,N));                 ...
       end                                                          
       for k=N-1:-1:1                                         
         for i=IstrU-1:Iend             
           P(i,j,k)=P(i,j,k+1)+                                        ...
                    HalfGRho*((rho(i,j,k+1)+rho(i,j,k)) .*               ...
                              (z_r(i,j,k+1)-z_r(i,j,k))-               ...
                              OneFifth*               ...
                              ((dR(i,k+1)-dR(i,k)) .*               ...
                               (z_r(i,j,k+1)-z_r(i,j,k)-               ...
                                OneTwelfth*                            ...
                                (dZ(i,k+1)+dZ(i,k)))-                  ...
                               (dZ(i,k+1)-dZ(i,k)) .*                  ...
                               (rho(i,j,k+1)-rho(i,j,k)-               ...
                                OneTwelfth*                            ...
                                (dR(i,k+1)+dR(i,k)))));
          end
        end
      end
%
%-----------------------------------------------------------------------
%%  Compute XI-component pressure gradient term.
%-----------------------------------------------------------------------
%
      for k=N:-1:1
        for j=Jstr:Jend
          for i=IstrU-1:Iend+1
            aux(i,j)=z_r(i,j,k)-z_r(i-1,j,k);
            FC(i,j)=rho(i,j,k)-rho(i-1,j,k);
          end
        end
%
        for j=Jstr:Jend
          for i=IstrU-1:Iend
            cff=2.0*aux(i,j)*aux(i+1,j);
            if (cff > eps) 
              cff1=1.0/(aux(i,j)+aux(i+1,j));
              dZx(i,j)=cff*cff1;
            else
              dZx(i,j)=0.0;
            end
            cff1=2.0*FC(i,j)*FC(i+1,j);
            if (cff1 > eps) 
              cff2=1.0/(FC(i,j)+FC(i+1,j));
              dRx(i,j)=cff1*cff2;
            else
              dRx(i,j)=0.0;
            end
          end
        end
%
        for j=Jstr:Jend
          for i=IstrU:Iend
            dPdx(i,j,k)=on_u(i,j)*0.5*                           ...
                          (Hz(i,j,k)+Hz(i-1,j,k))*                     ...
                          (P(i-1,j,k)-P(i,j,k)-                        ...
                           HalfGRho*                                   ...
                           ((rho(i,j,k)+rho(i-1,j,k))*                 ...
                            (z_r(i,j,k)-z_r(i-1,j,k))-                 ...
                             OneFifth*                                 ...
                             ((dRx(i,j)-dRx(i-1,j))*                   ...
                              (z_r(i,j,k)-z_r(i-1,j,k)-                ...
                               OneTwelfth*                             ...
                               (dZx(i,j)+dZx(i-1,j)))-                 ...
                              (dZx(i,j)-dZx(i-1,j))*                   ...
                              (rho(i,j,k)-rho(i-1,j,k)-                ...
                               OneTwelfth*                             ...
                               (dRx(i,j)+dRx(i-1,j))))));
          end
        end
      end
%
%-----------------------------------------------------------------------
%%  ETA-component pressure gradient term.
%-----------------------------------------------------------------------
%
      for k=N:-1:1
        for j=JstrV-1:Jend+1
          for i=Istr:Iend
            aux(i,j)=z_r(i,j,k)-z_r(i,j-1,k);
            FC(i,j)=rho(i,j,k)-rho(i,j-1,k);
          end
        end
%
        for j=JstrV-1:Jend
          for i=Istr:Iend
            cff=2.0*aux(i,j)*aux(i,j+1);
            if (cff > eps) 
              cff1=1.0/(aux(i,j)+aux(i,j+1));
              dZx(i,j)=cff*cff1;
            else
              dZx(i,j)=0.0;
            end
            cff1=2.0*FC(i,j)*FC(i,j+1);
            if (cff1 > eps) 
              cff2=1.0/(FC(i,j)+FC(i,j+1));
              dRx(i,j)=cff1*cff2;
            else
              dRx(i,j)=0.0;
            end
          end
        end
%
        for j=JstrV:Jend
          for i=Istr:Iend
            dPdy(i,j,k)=om_v(i,j)*0.5*                           ...
                          (Hz(i,j,k)+Hz(i,j-1,k))*                     ...
                          (P(i,j-1,k)-P(i,j,k)-                        ...
                           HalfGRho*                                   ...
                           ((rho(i,j,k)+rho(i,j-1,k))*                 ...
                            (z_r(i,j,k)-z_r(i,j-1,k))-                 ...
                             OneFifth*                                 ...
                             ((dRx(i,j)-dRx(i,j-1))*                   ...
                              (z_r(i,j,k)-z_r(i,j-1,k)-                ...
                               OneTwelfth*                             ...
                               (dZx(i,j)+dZx(i,j-1)))-                 ...
                              (dZx(i,j)-dZx(i,j-1))*                   ...
                              (rho(i,j,k)-rho(i,j-1,k)-                ...
                               OneTwelfth*                             ...
                               (dRx(i,j)+dRx(i,j-1))))));
          end
        end
      end
 OneTwelfth
