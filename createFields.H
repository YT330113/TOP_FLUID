/***********************************************************/
/**********************创建设计变量场**********************/
/***********************************************************/

Info << "Reading optimization parameters\n" << endl;

/********************从optProperties字典文件中读取参数*****************/

IOdictionary optProperties                      //创建字典，用于输入输出
(
   IOobject
   (
       "optProperties",                                     // 字典名称
       runTime.constant(),                              // 字典路径
       mesh,
       IOobject::MUST_READ_IF_MODIFIED, //字典读控制（从字典文件读取到程序）
       IOobject::NO_WRITE                                   // 字典写控制（程序写入到字典文件）
   )
);
scalar voluse(readScalar(optProperties.lookup("voluse")));      // 体积约束
scalar GeoDim(readScalar(optProperties.lookup("GeoDim"))); // 模型维度，三维/二维
scalar PowerDiss0(readScalar(optProperties.lookup("PowerDiss0"))); // 能量耗散归一化值

Switch solid_area(optProperties.lookup("solid_area"));              //bool型变量，判断是否存在不可设计域
Switch fluid_area(optProperties.lookup("fluid_area"));              // 流体的不可设计域
Switch test_area(optProperties.lookup("test_area"));                    //用于测试能量耗散

scalar coff_PDE(readScalar(optProperties.lookup("coff_PDE")));//用于PDE过滤，coff_PDE * len=r,r为coff_PDE倍的单元长度;coff_PDE可取1.0-2.0
dimensionedScalar b("b", dimensionSet(0, -2, 0, 0, 0, 0, 0), 1);// b就是PDE求导数中的1/( r^2 ),后续会定义具体值  

//-------------------------------------------------------------------------------//

/**************************定义伪密度场gamma**************************/

// xp是PDE过滤的；xh是Heviside过滤的，xh是最终更新的设计变量，
//用于下一次迭代时控制方程的求解 ; x是优化变量，即MMA优化求解器中的当前点
Info << "Reading pseudo-density field\n"<< endl;
volScalarField xp                           //创建体标量场
(
    IOobject
    (
        "xp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    scalar(voluse),
    // 如果在算例0文件夹下存在gamma的初边值文件，则读入文件中的值
    //作为设计变量的初值；否则采用voluse作为初值。
    zeroGradientFvPatchScalarField::typeName   
    // 指定场边界信息，所有边界均为zeroGradient 边界
    // 体心标量场类型的数据不仅包含标量场在控制体积中心存储的离散值，
    // 还包含了该场的边界信息
);
volScalarField xh
(
    IOobject
    (
        "xh",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    scalar(voluse),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField x(xh);
volScalarField drho(xh);                            //这个drho是xh对xp的导数
//--------------------------------------------------------------------------------------------//


/*****************************定义体积约束灵敏度****************************/

volScalarField gsensVol         //约束函数灵敏度，对x的导数                              
(
    IOobject
    (
        "gsensVol",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    xh,
    zeroGradientFvPatchScalarField::typeName
);
volScalarField gsenshVol(gsensVol); //约束函数灵敏度，对xh的导数
//--------------------------------------------------------------------------------------------//


/*****************************定义非设计域****************************/

labelList cells_solid, cells_fluid, cells_test; //定义标签列表，存放非设计域单元编号

if (solid_area)
{
  word zoneName="solid_area";
  label zoneID=mesh.cellZones().findZoneID(zoneName); //找到solid_area cellzone的序号，
                                                                            //注意是整个solid_area的序号，不是每个单元的序号
  cells_solid=mesh.cellZones()[zoneID];       //将非设计域单元编号存在cells列表里

  forAll(cells_solid, i)
  {
     x[cells_solid[i]]=0;                                           //对非设计域单元将设计变量置1
  }
}

if(fluid_area)
{
  word zoneName="fluid_area";
  label zoneID=mesh.cellZones().findZoneID(zoneName);
  cells_fluid=mesh.cellZones()[zoneID];
  forAll(cells_fluid, i)
  {
     x[cells_fluid[i]]=1.0;                                                          // 给所有fluid单元伪密度x赋值为1
  }
}

  if(test_area)                 // test单元用于测试能量耗散
{
  word zoneName="test_area";
  label zoneID=mesh.cellZones().findZoneID(zoneName);
  cells_test=mesh.cellZones()[zoneID];
  forAll(cells_test, i)
  {
     x[cells_test[i]]=1.0;                                                         
     xh[cells_test[i]]=1.0;
  }
}
//---------------------------------------------------------------------------------------------//


/* **************************查找边界处的网格用于计算定义在边界上的函数值 **************************/ 

dictionary conFunc = mesh.solutionDict().subDict("constraintFunctionDict_flow");  // 查找fvSolution字典下的“constraintFunctionDict_flow”栏
int nObjPatch =conFunc.lookupOrDefault<scalar>("numberConstraintPatches",2);   // 指目标边界面的个数，默认就是inlet和outlet两个面
wordList conPatchNames = conFunc.lookup("constraintPatchesNames");   // 边界面的名字,conPatchNames是个列表，其值为(inlet outlet等等)
label conPatchList[nObjPatch];  //一个label数组，有两个:conPatchList[0],conPatchList[1]，conPatchList在“costFunction”计算能量耗散用到
int iLoop;
for(iLoop=0; iLoop<nObjPatch; iLoop++)  // 这一步是为了获得需要计算能力那个耗散的patch的ID号
{
    conPatchList[iLoop]=mesh.boundaryMesh().findPatchID(conPatchNames[iLoop]);// mesh.boundaryMesh().findPatchID是查找某个大面的ID，
                                                                                                                                                                            //比如说inlet的ID是0，outlet的ID是1...
}  
//------------------------------------------------------------------------------------------------------------------------------------------//
