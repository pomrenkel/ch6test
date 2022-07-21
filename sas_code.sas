SAS/IML Code to Generate Simulation Data

options nodate nonumber linesize=90;
proc iml;
goseed1=3001111;
reps=3000;
goal=1000;
factors=3; p=6; vars=factors*p;
half=p/2;
load={0.7,0.4};
error={0.51,0.84};
phi={1.0 0.4 0.4,
     0.4 1.0 0.4,
	 0.4 0.4 1.0};
n={100,200,500,1000};
keepthis=repeat(0,nrow(n),8);
storthis=repeat(0,goal*nrow(n),10);
lambda=repeat(0,vars,factors);
thetdelt=repeat(0,vars,vars);
sds=repeat(0,vars,vars);

do f=1 to factors;
  do h=1 to half;
    lambda[p*(f-1)+h,f]=load[1,1];
    lambda[p*(f-1)+h+half,f]=load[2,1];
    thetdelt[p*(f-1)+h,p*(f-1)+h]=error[1,1];
	thetdelt[p*(f-1)+h+half,p*(f-1)+h+half]=error[2,1];
  end;
  do j=1 to p;
    sds[p*(f-1)+j,p*(f-1)+j]=root(lambda[p*(f-1)+j,f]*lambda[p*(f-1)+j,f]*phi[f,f]+thetdelt[p*(f-1)+j,p*(f-1)+j]);
  end;
end;

sigma=lambda*phi*lambda`+thetdelt;
ds=det(sigma);
R1=inv(sds)*sigma*inv(sds);
g=root(R1);

do sampsize=1 to nrow(n); /** Begin sample size loop **/

seed1=goseed1+20000*(sampsize-1);
holdthis=repeat(0,goal,8);
w=0; added=0;
do i=1 to reps; /** Begin replication loop **/
z1=normal(repeat(seed1,n[sampsize,1],vars));
D=z1*g*sds;

file 'C:\research\chapter\runthis.dat';
  do r=1 to n[sampsize,1];
    do c=1 to vars;
      put (D[r,c]) +1 @;
	end;
    put;
  end;
closefile 'C:\research\chapter\runthis.dat';

run system('C:\research\chapter\goforth3');
infile 'C:\research\chapter\pf06_mot.out';
input / / / / / / / / / / / / / / / / / / / / / / / / /
      / / / / / / / / / / / / / / / / / / / / / / / / /
      / / / / / / / / / / / @7 test $char2.; *print test;
if test="OF" then do;
input / / / @42 chisq;
input / @45 chi_pval;
input / / / / / / / / / @46 CFI;
input @46 TLI;
input / / / / / / / / / / / / / / / / @46 RMSEA;
input @46 RMSEA_LB @53 RMSEA_UB;
input / / / / @46 SRMR;
w=w+1;
holdthis[w,1]=chisq; holdthis[w,2]=chi_pval; holdthis[w,3]=CFI; holdthis[w,4]=TLI;
holdthis[w,5]=RMSEA; holdthis[w,6]=RMSEA_LB; holdthis[w,7]=RMSEA_UB; holdthis[w,8]=SRMR;
do j=1 to 8;
  storthis[goal*(sampsize-1)+w,1]=n[sampsize,1];
  storthis[goal*(sampsize-1)+w,j+1]=holdthis[w,j];
  storthis[goal*(sampsize-1)+w,10]=added;
end;
end; /* If worked do loop */
else added=added+1;
closefile 'C:\research\chapter\pf06_mot.out';
if w=goal then i=reps;

end; /** Replication loop **/

do j=1 to 8;
  keepthis[sampsize,j]=holdthis[:,j];
end;

end; /** Sample size loop **/

file 'C:\research\chapter\keepthis.dat';
  do r=1 to nrow(n);
    do c=1 to 8;
      put (keepthis[r,c]) +1 @;
	end;
    put;
  end;
closefile 'C:\research\chapter\keepthis.dat';
file 'C:\research\chapter\storthis.dat';
  do r=1 to nrow(n)*goal;
    do c=1 to ncol(storthis);
      put (storthis[r,c]) +1 @;
	end;
    put;
  end;
closefile 'C:\research\chapter\storthis.dat';

quit;
