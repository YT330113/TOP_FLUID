/*************************************************************************************************/
/***********************灵敏度分析及灵敏度链式推导*********************************************/
/*************************************************************************************************/


Info<<"\nsensitivity analysis\n"<<endl;

/**************************计算柔度目标函数和体积约束的灵敏度****************************/

fsenshPowerDiss=-alphaMax*(1+qu)*qu/((qu+xh)*(qu+xh))*(U&Ua);
for(i=0;i<n;i++)
{
 fsenshPowerDiss[i]=fsenshPowerDiss[i]*mesh.V()[i];
}

//灵敏度归一化
fsenshPowerDiss=fsenshPowerDiss/gMax(mag(fsenshPowerDiss.primitiveField()));//mag()返回绝对值,// gMax返回内部和边界场的最大值
gsenshVol.primitiveFieldRef()=mesh.V()/gMax(mesh.V()); // primitiveFieldRef()返回无量纲的 内部场，
//mag()取模，上面两行返回了各自灵敏度初以绝对值最大的那个，为什么？？？归一化
//-------------------------------------------------------------------------------------------------------------------------//

/****************************************不可设计域的灵敏度置零********************************/

if(solid_area)
{
   forAll(cells_solid, i)
   {
      fsenshPowerDiss[cells_solid[i]]=0;
      gsenshVol[cells_solid[i]]=0;
   }
}
if(fluid_area)
{
   forAll(cells_fluid, i)
   {
      fsenshPowerDiss[cells_fluid[i]]=0;
      gsenshVol[cells_fluid[i]]=0;
   }
}
//------------------------------------------------------------------------------------------------------------------------------// 

/***********************************链式求对设计变量x的导数用于MMA****************************/

//这是在计算sp的表达式，函数phi对xh的导数乘以xh对xp的导数

/*if(opt>1)
{
  for(i=0;i<n;i++)
  {
     if(xp[i]<=eta5)
     {
         drho[i] = del * Foam::exp(-del * (1 - xp[i] / eta5)) + Foam::exp(-del);
         // del就是heviside过滤里的惩罚强度参数betta
         // eta5就是截断值小于该值的伪密度会倾向于罚为0；大于该值的倾向于罚为1。
         //每个优化步都需要计算 eta值来保证过滤前后材料用量一致，即保体积性
     }
     else
     {
         drho[i] = del * Foam::exp(-del * (xp[i] - eta5) / (1 - eta5)) + Foam::exp(-del);
     }
  }
  fsenshPowerDiss = fsenshPowerDiss * drho; // sp的表达式，函数phi对xh的导数乘以xh对xp的导数
  gsenshVol = gsenshVol * drho;
}




// 求解PDE得到函数phi对设计变量x的导数fens，这个导数用于MMA求解器

solve(fvm::laplacian(fsensPowerDiss) - fvm::Sp(b, fsensPowerDiss) + fsenshPowerDiss * b);//b就是1/r^2，之前已经定义
solve(fvm::laplacian(gsensVol) - fvm::Sp(b, gsensVol) + gsenshVol * b);
*/
//不过滤
fsensPowerDiss = fsenshPowerDiss ; // sp的表达式，函数phi对xh的导数乘以xh对xp的导数
gsensVol = gsenshVol;
//最终得到了函数phi对设计变量x的导数，可带入MMA优化求解器求解更新设计变量
//-------------------------------------------------------------------------------------------------------------------------//