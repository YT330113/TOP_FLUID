    //**************************************动量预测***************************************//
//for(i=0;i<2500;i++)//执行40次流场计算循环，ajointShapeOptimizationFoam没有循环
{
    //Info << "\nsimple Iter: " << i << nl;

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)                                                            // 对流项目
    + turbulence->divDevReff(U)                                   // 扩散项
    + fvm::Sp(alpha, U)                                                       // 人工阻力源项
    ==
        fvOptions(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();               //对 UEqn的系数矩阵进行一些调整，以使它更满足对角占优。 fvSolution 里面如果没有定义松弛因子，则默认值为1。

    fvOptions.constrain(UEqn);
    solve(UEqn == -fvc::grad(p));
    fvOptions.correct(U);
//*************************************************压力修正*****************************************//
volScalarField rAU(1.0/UEqn.A());//计算1/aP，其中UEqn.A()获得Ap，从矩阵的角度看即是由系数矩阵的对角线元素组成的矩阵
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
tUEqn.clear();
surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));//计算HbyA，即式(2.9)至(2.10)的过程
adjustPhi(phiHbyA, U, p);               //对流率进行调整，以保证方程的连续性

tmp<volScalarField> rAtU(rAU);//可能与simpleC有关

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, U, phiHbyA, rAU);//更新压力边界条件以保证通量计算守恒

// Non-orthogonal pressure corrector loop
while (simple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAU, p) == fvc::div(phiHbyA)//组建压力泊松方程
    );

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve();//求解压力泊松方程

    if (simple.finalNonOrthogonalIter())
    {
        phi = phiHbyA - pEqn.flux();          // fvc::flux(U)  // 通量的值， 该句等同于fvc::interpolate(U) & mesh.Sf()
        //在最后一次非正交修正结束后更新速度，即计算式(2.15)或式(2.24)，但这里更新的是速度的通量，而不是速度本身
    }
}

#include "continuityErrs.H"//计算连续性方程计算误差

// Explicitly relax pressure for momentum corrector
p.relax();//压力场进行显示亚松弛

// Momentum corrector动量修正
U = HbyA - rAU*fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);

//湍流修正
laminarTransport.correct();
turbulence->correct();

//runTime.write();

}



