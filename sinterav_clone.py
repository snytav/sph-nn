import numpy as np
import pandas as pd

d = {'GamE': [978030.00000, 978049.00000,978031.80000,], 'col2': [3, 4]}


NormGrav: Array[0.
.4] of
TNormGrav =
((GamE: 978030.00000; Bt: 0.0053020; Bt1: 0.0000070; C20: -1083.46E-06; C40: 2.72E-06; aE: 0),
 (GamE: 978049.00000; Bt: 0.0052884; Bt1: 0.0000059; C20: -1091.87E-06; C40: 2.42E-06; aE: 6378172),
 (GamE: 978031.80000; Bt: 0.0053024; Bt1: 0.0000059; C20: -1082.78E-06; C40: 2.37E-06; aE: 6378172),
 (GamE: 978031.85000; Bt: 0.0053024; Bt1: 0.0000059; C20: -1083.63E-06; C40: 1.62E-06; aE: 6378172),
 (GamE
  : 978032.53359; Bt:         0; Bt1:         0; C20: -1082.63E-06; C40: 2.37E-06; aE: 6371008.7714; e: 8.1819190842622E-2));


def fLGN1(M, N, X):
# Var
#      A,S,U,F,D,P: Extended;
# begin
   K_LGN1 = 0
   H_LGN1 = 0
   B_LGN1 = 0
   C_LGN1 = 0
   if N == 0:
      C_LGN1 = np.sqrt(0.5)
      K_LGN1 = 0
      B_LGN1 = 0.0
      H_LGN1 = C_LGN1
      LGN1   = C_LGN1

   else:
     if M > K_LGN1:
          C_LGN1 = (2*N+1)*(1-X*X)/(2*M)
          C_LGN1 = np.sqrt(C_LGN1)*H_LGN1
          K_LGN1 = M
          B_LGN1 = 0
          H_LGN1 = C_LGN1
          LGN1   = C_LGN1
     else:
          A = 2*N
          S = np.sqrt(A*A-1)
          U = N*N-M*M
          F = 1.0/np.sqrt(U)
          P = np.sqrt(np.abs((A+1)/(A-3)*((N-1)*(N-1)-M*M)))
          D = (X*C_LGN1*S-B_LGN1*P)*F
          B_LGN1 = C_LGN1
          C_LGN1 = D
          LGN1   = C_LGN1
   return LGN1

def NormGamm(GT,S2Fi, S22Fi):
    return NormGrav[T].GamE*(1+NormGrav[GT].Bt*S2Fi-NormGrav[GT].Bt1*S22Fi)


def VBV4(MR,VFlag,GamT,FI,AL,N0,K_vbv,B_vbv):
# var
#     CA, SA, P1, R0, R1, SinFi,SQRT2,SQRT05: Extended;
# 	  m, n: Integer;
#     B_LGN1, C_LGN1, H_LGN1: Extended;
#     K_LGN1: Integer;
#     RadiusGamm: Extended;
#     GMFi,S2Fi,S22Fi: Extended;

# var
#   tfOut: TextFile;
#   c_FNAME,s1,s2: String;
   RadiusGamm = 1
   if VFlag:
      RadiusGamm = MR
   sinFI  = np.sin(FI)
   SQRT2  = np.sqrt(2)
   SQRT05 = np.sqrt(0.5)
   S2Fi   = sinFI*sinFI
   S22Fi  = np.sin(2*FI)*np.sin(2*FI)
   GMFi   = NormGamm(GamT,S2Fi,S22Fi)

   if GamT == 4:
      GMFi  = NormGrav[GamT].GamE*(1+0.001931852654*np.sin(Fi)*np.sin(Fi))/np.sqrt(1-0.0066943799*np.sin(fi)*np.sin(fi))
      fi    = np.arctan((np.tan(fi)*(1-(1.0/298.257))*(1.0-(1.0/298.257))))
      sinFI = np.sin(FI)


     K_vbv[2,0]:= BaseC20-NormGrav[GamT].C20/SQRT(5);
     K_vbv[4,0]:= BaseC40-NormGrav[GamT].C40/SQRT(9);

   // Set the name of the file that will be created
   s1 := floattostr(fi);
   s2 := floattostr(al);
   s1 := concat(s1,'_');
   s2 := concat(s2,'_.txt');
   c_fname:= concat('sphar_',concat(s1,s2));

   AssignFile(tfOut, C_FNAME);
   rewrite(tfOut);


   for m := 0 to N0 do
    begin
     R0 := 0;
     R1 := 0;
     CA := cos(M*AL);
     SA := sin(M*AL);
     for n:= m to N0 do
      begin
       P1:= LGN1(m,n, sinFI);
       If m = 0 Then p1:= p1*SQRT2
                Else p1:= p1*2;
       If n<2 Then Continue;
       If VFlag Then P1:= P1*RadiusGamm
                Else P1:= P1*(GMFi*(n-1));
       If m <> 0 Then R1:= R1+K_vbv[m-1,n]*P1;
                      R0:= R0+K_vbv[n,m]*P1;

       writeln(tfOut,'m,n,ca,r0,sa,r1,K_vbv[n,m],P1 ',
                      m:5,n:5,ca:15:5,r0:15:5,sa:15:5,r1:15:5,K_vbv[n,m]:15:5,P1:25:10);

      end;
     B_vbv := B_vbv + CA*R0+SA*R1;
       end;
    CloseFile(tfOut);


  K_vbv[2,0]:= BaseC20;
  K_vbv[4,0]:= BaseC40;
end;
