/***********************Adjoint Pressure-velocity SIMPLE corrector************************/

//for(i=0;i<2500;i++)//执行30次伴随场计算循环，ajointShapeOptimizationFoam没有循环
{
    //Info << "\nsimple Iter: " << i << nl;
    // Adjoint Momentum predictor
    volVectorField adjointTransposeConvection((fvc::grad(Ua) & U));

    tmp<fvVectorMatrix> tUaEqn
    (
        fvm::div(-phi, Ua)
        - adjointTransposeConvection
        + turbulence->divDevReff(Ua)
        + fvm::Sp(alpha, Ua)
        ==
        fvOptions(Ua)
    );
    fvVectorMatrix& UaEqn = tUaEqn.ref();

    UaEqn.relax();

    fvOptions.constrain(UaEqn);

    solve(UaEqn == -fvc::grad(pa));

    fvOptions.correct(Ua);
//****************************************
    volScalarField rAUa(1.0/UaEqn.A());
    volVectorField HbyAa("HbyAa", Ua);
    HbyAa = rAUa*UaEqn.H();
    tUaEqn.clear();
    surfaceScalarField phiHbyAa("phiHbyAa", fvc::flux(HbyAa));
    adjustPhi(phiHbyAa, Ua, pa);

    // Update the pressure BCs to ensure flux consistency
    //constrainPressure(pa, Ua, phiHbyAa, rAUa);   //这句原本没有

    // Non-orthogonal pressure corrector loop
    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix paEqn
        (
            fvm::laplacian(rAUa, pa) == fvc::div(phiHbyAa)
        );

        paEqn.setReference(paRefCell, paRefValue);
        paEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            phia = phiHbyAa - paEqn.flux();
        }
    }

    #include "adjointContinuityErrs.H"

    // Explicitly relax pressure for adjoint momentum corrector
    pa.relax();

    // Adjoint momentum corrector
    Ua = HbyAa - rAUa*fvc::grad(pa);
    Ua.correctBoundaryConditions();
    fvOptions.correct(Ua);

        //Ua.storePrevIter();         //不知道干嘛的
        //pa.storePrevIter();     
        //phia.storePrevIter();
}
