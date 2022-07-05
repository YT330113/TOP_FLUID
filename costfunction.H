/**************************************************************************/
/*****************目标函数和约束函数值计算*******************************/
/**************************************************************************/


Info <<"\nCaculate cost function\n"<<endl;
/****************************计算能量耗散和体积约束******************************/

PowerDiss=0;
for (i=0; i<nObjPatch; i++)
{
     PowerDiss =PowerDiss - sum(phi.boundaryField()[conPatchList[i]]*(p.boundaryField()[conPatchList[i]]+ 0.5*magSqr(U.boundaryField()[conPatchList[i]]))); 
}

V = fvc::domainIntegrate(xh).value() / area - voluse; // g(x)<=0;     area和voluse都是无量刚标量
//-----------------------------------------------------------------------------------------------------------//

/****************************输出每一迭代步的目标函数值*********************/

if (Pstream::master())
{
    ofstream outfile1("Voluse.txt", std::ios::app);
    outfile1 <<runTime.timeName()<<"\t"<< V << "\n";
    outfile1.close();

    ofstream outfile2("PowerDissRatio.txt",std::ios::app);
    outfile2<<runTime.timeName()<<"\t"<<PowerDiss<<"\n";
    outfile2.close(); 

    ofstream outfile3("Time.txt", std::ios::app);
    outfile3 <<runTime.timeName()<<"\t"<< runTime.clockTimeIncrement() << "\n";
    outfile3.close();    
}
