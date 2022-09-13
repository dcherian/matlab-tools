function [gradPx,gradPy] = roms_prsgd31(K)

%---------------------------------------------------------------------------
% Compute pressure gradient components.
%---------------------------------------------------------------------------
%
% NOTE: fac2=0 because the balanced component should consist of the
% baroclinic pressure gradient only.

fac1=0.5*K.g/K.rho0;
fac2=0.0;                        % originally, fac2=g (not used anyway)
fac3=0.25*K.g/K.rho0;

% Compute surface baroclinic pressure gradient (phix, m2/s2) and
% its vertically integrated values (phix_bar, m3/s2; gradPx, m/s)
% in the XI-direction (at U-points).
%
%   (i-1,j,N)      1:L ,1:M,N         DO j=Jstr-1,Jend
%   (i  ,j,N)      2:Lp,1:M,N           DO i=Istr,Iend+1 

cff1(1:L,1:M)=K.Zw(2:Lp,1:M,Np)-                                         ...
              K.Zr(2:Lp,1:M,N )+                                         ...
              K.Zw(1:L ,1:M,Np)-                                         ...
              K.Zr(1:L ,1:M,N );

phix(1:L,1:M)=fac1.*(deltaR_b(2:Lp,1:M,N)-                               ...
                     deltaR_b(1:L ,1:M,N)).*cff1(1:L,1:M);

phix_bar(1:L,1:M)=0.5.*(K.Hz(1:L ,1:M,N)+K.Hz(2:Lp,1:M,N)).*             ...
                  phix(1:L,1:M).*                                        ...
                  K.umask(1:L,1:M);

gradPx(1:L,1:M,N)=0.5.*phix(1:L,1:M).*                                   ...
                       (K.pm(1:L ,1:M)+K.pm(2:Lp,1:M))./                 ...
                       ( K.f(1:L ,1:M)+ K.f(2:Lp,1:M));

% Compute interior baroclinic pressure gradient (phix, m2/s2) and
% its vertically integrated values (phix_bar, m3/s2; gradPx, m/s)
% in the XI-direction (at U-points). Differentiate and then
% vertically integrate.
%
%   (i-1,j,k  )    1:L ,1:M,k        DO k=1,N-1
%   (i  ,j,k  )    2:Lp,1:M,k          DO j=Jstr-1,Jend
%   (i-1,j,k+1)    1:L ,1:M,k+1          DO i=Istr,Iend+1
%   (i  ,j,k+1)    2:Lp,1:M,k+1

for k=Nm:-1:1,

  cff1(1:L,1:M)=1.0./((K.Zr(2:Lp,1:M,k+1)-K.Zr(2:Lp,1:M,k  )).*          ...
                      (K.Zr(1:L ,1:M,k+1)-K.Zr(1:L ,1:M,k  )));

  cff2(1:L,1:M)=K.Zr(2:Lp,1:M,k  )-K.Zr(1:L ,1:M,k  )+                   ...
                K.Zr(2:Lp,1:M,k+1)-K.Zr(1:L ,1:M,k+1);

  cff3(1:L,1:M)=K.Zr(2:Lp,1:M,k+1)-K.Zr(2:Lp,1:M,k  )-                   ...
                K.Zr(1:L ,1:M,k+1)+K.Zr(1:L ,1:M,k  );

  gamma(1:L,1:M)=0.125.*cff1.*cff2.*cff3;

  cff1(1:L,1:M)=(1.0+gamma(1:L,1:M)).*(deltaR_b(2:Lp,1:M,k+1)-           ...
                                       deltaR_b(1:L ,1:M,k+1))+          ...
                (1.0-gamma(1:L,1:M)).*(deltaR_b(2:Lp,1:M,k  )-           ...
                                       deltaR_b(1:L ,1:M,k  ));

  cff2(1:L,1:M)=deltaR_b(2:Lp,1:M,k+1)+deltaR_b(1:L ,1:M,k+1)-           ...
                deltaR_b(2:Lp,1:M,k  )-deltaR_b(1:L ,1:M,k  );

  cff3(1:L,1:M)=K.Zr(2:Lp,1:M,k+1)+K.Zr(1:L ,1:M,k+1)-                   ...
                K.Zr(2:Lp,1:M,k  )-K.Zr(1:L ,1:M,k  );

  cff4(1:L,1:M)=(1.0+gamma(1:L,1:M)).*(K.Zr(2:Lp,1:M,k+1)-               ...
                                       K.Zr(1:L ,1:M,k+1))+              ...
                (1.0-gamma(1:L,1:M)).*(K.Zr(2:Lp,1:M,k  )-               ...
                                       K.Zr(1:L ,1:M,k  ));

  phix(1:L,1:M)=phix(1:L,1:M)+fac3.*(cff1.*cff3-cff2.*cff4);

  phix_bar(1:L,1:M)=phix_bar(1:L,1:M)+                                   ...
                    0.5.*(K.Hz(1:L ,1:M,k)+K.Hz(2:Lp,1:M,k)).*           ...
                    phix(1:L,1:M).*                                      ...
                    K.umask(1:L,1:M);

  gradPx(1:L,1:M,k)=0.5.*phix(1:L,1:M).*                                 ...
                         (K.pm(1:L ,1:M)+K.pm(2:Lp,1:M))./               ...
                         ( K.f(1:L ,1:M)+ K.f(2:Lp,1:M));

end,

% Compute surface baroclinic pressure gradient (phie, m2/s2) and
% its vertically integrated value (phie_bar, m3/s2; gradPy, m/s)
% in the ETA-direction (at V-points).
%
%   (i,j-1,N)      1:L,1:M ,N         DO j=Jstr,Jend+1
%   (i,j  ,N)      1:L,2:Mp,N           DO i=Istr-1,Iend

cff1(1:L,1:M)=K.Zw(1:L,2:Mp,Np)-                                         ...
              K.Zr(1:L,2:Mp,N )+                                         ...
              K.Zw(1:L,1:M ,Np)-                                         ...
              K.Zr(1:L,1:M ,N );
  
phie(1:L,1:M)=fac1.*(deltaR_b(1:L,2:Mp,N)-                               ...
                     deltaR_b(1:L,1:M ,N)).*cff1(1:L,1:M);

phie_bar(1:L,1:M)=0.5.*(K.Hz(1:L,1:M ,N)+K.Hz(1:L,2:Mp,N)).*             ...
                  phie(1:L,1:M).*                                        ...
                  K.vmask(1:L,1:M);

gradPy(1:L,1:M,N)=0.5.*phie(1:L,1:M).*                                   ...
                       (K.pn(1:L,1:M )+K.pn(1:L,2:Mp))./                 ...
                       ( K.f(1:L,1:M )+ K.f(1:L,2:Mp));

% Compute interior baroclinic pressure gradient (phie, m2/s2) and
% its vertically integrated value (phie_bar, m3/s2; gradPy, m/s)
% in the ETA-direction (at V-points). Differentiate and then
% vertically integrate.
%
%   (i,j-1,k  )    1:L,1:M ,k         DO k=1,N-1
%   (i,j  ,k  )    1:L,2:Mp,k           DO j=Jstr,Jend+1
%   (i,j-1,k+1)    1:L,1:M ,k+1           DO i=Istr-1,Iend
%   (i,j  ,k+1)    1:L,2:Mp,k+1

for k=Nm:-1:1,

  cff1(1:L,1:M)=1.0./((K.Zr(1:L,2:Mp,k+1)-K.Zr(1:L,2:Mp,k  )).*          ...
                      (K.Zr(1:L,1:M ,k+1)-K.Zr(1:L,1:M ,k  )));
    
  cff2(1:L,1:M)=K.Zr(1:L,2:Mp,k  )-K.Zr(1:L,1:M ,k  )+                   ...
                K.Zr(1:L,2:Mp,k+1)-K.Zr(1:L,1:M ,k+1);

  cff3(1:L,1:M)=K.Zr(1:L,2:Mp,k+1)-K.Zr(1:L,2:Mp,k  )-                   ...
                K.Zr(1:L,1:M ,k+1)+K.Zr(1:L,1:M ,k  );
   
  gamma(1:L,1:M)=0.125.*cff1.*cff2.*cff3;

  cff1(1:L,1:M)=(1.0+gamma(1:L,1:M)).*(deltaR_b(1:L,2:Mp,k+1)-           ...
                                       deltaR_b(1:L,1:M ,k+1))+          ...
                (1.0-gamma(1:L,1:M)).*(deltaR_b(1:L,2:Mp,k  )-           ...
                                       deltaR_b(1:L,1:M ,k  ));

  cff2(1:L,1:M)=deltaR_b(1:L,2:Mp,k+1)+deltaR_b(1:L,1:M ,k+1)-           ...
                deltaR_b(1:L,2:Mp,k  )-deltaR_b(1:L,1:M ,k  );

  cff3(1:L,1:M)=K.Zr(1:L,2:Mp,k+1)+K.Zr(1:L,1:M ,k+1)-                   ...
                K.Zr(1:L,2:Mp,k  )-K.Zr(1:L,1:M ,k  );

  cff4(1:L,1:M)=(1.0+gamma(1:L,1:M)).*(K.Zr(1:L,2:Mp,k+1)-               ...
                                       K.Zr(1:L,1:M ,k+1))+              ...
                (1.0-gamma(1:L,1:M)).*(K.Zr(1:L,2:Mp,k  )-               ...
                                       K.Zr(1:L,1:M ,k  ));

  phie(1:L,1:M)=phie(1:L,1:M)+                                           ...
                fac3.*(cff1.*cff3-cff2.*cff4);

  phie_bar(1:L,1:M)=phie_bar(1:L,1:M)+                                   ...
                    0.5.*(K.Hz(1:L,1:M ,k)+K.Hz(1:L,2:Mp,k)).*           ...
                    phie(1:L,1:M).*                                      ...
                    K.vmask(1:L,1:M);

  gradPy(1:L,1:M,k)=0.5.*phie(1:L,1:M).*                                 ...
                         (K.pn(1:L,1:M )+K.pn(1:L,2:Mp))./               ...
                         ( K.f(1:L,1:M )+ K.f(1:L,2:Mp));

end,

clear cff1 cff2 cff3 cff4 gamma phie phix