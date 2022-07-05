/*************************************************************************************************/
/***************************************拓扑变量求解及过滤**************************************/
/*************************************************************************************************/


/********************************MMA求解器拓扑变量求解***********************************/

std::cout << "\nRun MMA algorithm\n"<< std::endl;

g[0] = V; // constraint function
for (i = 0; i < n; i++)
{
	df[i] = fsensPowerDiss[i] / n; // sensitivities of objective function
	dg[i] = gsensVol[i] / n; // sensitivities of constraint function
	}  

	// Set outer move limits
	for (int i=0;i<n;i++) {
		xmax[i] = Min(Xmax, xmma[i] + movlim);
		xmin[i] = Max(Xmin, xmma[i] - movlim);
	}

	// Call the update method
	mma->Update(xmma,df,g,dg,xmin,xmax);

	// Compute infnorm on design change
	double ch = 0.0;
	for (int i=0;i<n;i++) {
		ch = Max(ch,Abs(xmma[i]-xold[i]));
		xold[i] = xmma[i];
		x[i]=xmma[i];
	}

	// Print to screen
	std::cout << "\n设计变量的最大变化量:  "<< ch << std::endl;

	// Print to file
  if (Pstream::master())
{
    ofstream outfile4("Change.txt", std::ios::app);
    outfile4 << ch << "\n";
    outfile4.close();
}
//-------------------------------------------------------------------------------------------------------------------------//

//没有过滤
xh = x;

/********************************PDE过滤及Heviside过滤***********************************/

/***********************************PDE filter**********************************************/
/*
Info<<"\nPDE过滤\n"<<endl;
solve(fvm::laplacian(xp) - fvm::Sp(b, xp) + x * b);*/
//-------------------------------------------------------------------------------------------------------------------------//

//***********************************Heaviside_filter***************************************//
/*
Info<<"\nHeaviside过滤\n"<<endl;
del = Foam::min(0.2 * opt, 100);
eta0 = 0.0001;
eta1 = 1.0;
y0 = diff(xp, mesh.V(), del, eta0, n);
reduce(y0, sumOp<scalar>());

do
{
  eta5=(eta0+eta1)/2.0;
  y5=diff(xp,mesh.V(),del,eta5,n);
  reduce(y5, sumOp<scalar>());
  if(y0*y5<0)
  { 
     eta1=eta5;
  }
  else
  { 
     eta0=eta5;
     y0=y5;
  }
} while ((eta1-eta0)>0.0001);

for(i=0;i<n;i++)
{
   if(xp[i]<=eta5)
   {
     xh[i]=eta5*(Foam::exp(-del*(1-xp[i]/eta5))-(1-xp[i]/eta5)*Foam::exp(-del)); 
   }
   else
   {
     xh[i]=eta5+(1-eta5)*(1-Foam::exp(-del*(xp[i]-eta5)/(1-eta5))+(xp[i]-eta5)*Foam::exp(-del)/(1-eta5));
   }   
}
*/
//-------------------------------------------------------------------------------------------------------------------------//




