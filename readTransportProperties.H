/* *********************************************读取流体属性*********************************************** */
Info<< "\nReading TransportProperties\n" << endl;
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);
dimensionedScalar nu(transportProperties.lookup("nu"));                                     //流体粘度
dimensionedScalar alphaMax(transportProperties.lookup("alphaMax"));        //初始流阻系数的最大值 （刚开始优化时放松），alphaMax=2.5e3
dimensionedScalar alphaMAX(transportProperties.lookup("alphaMAX"));      // 流阻系数的最大值alpha_max,alphaMAX=2.5e7
dimensionedScalar alphamax(transportProperties.lookup("alphaMax"));        // 用于放松，alphaMax=alphamax*(opt/7.0+1.0);
scalar qu(0.01);                                                                                                      // 定义流阻系数插值常数q为0.01
volScalarField alpha(alphaMax*qu*(1-xh)/(qu+xh));                                                      // 体标量场流阻系数插值alpha_f

/**********************************************************************************************************************/
#include "createMRF.H"                                                      // simpleFoam有
Info<< "\nReading field p\n" << endl;
volScalarField p                                                                       // 定义压力场p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
Info<< "\nReading field U\n" << endl;                                 // 定义速度场U
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
#include "createPhi.H"                                                          // 定义通量场phi，面心标量场，48行到57是simpleFoam有的
label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, simple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());
singlePhaseTransportModel laminarTransport(U, phi);
autoPtr<incompressible::turbulenceModel> turbulence// turbulence是个autoPtr，turbulence与粘度设置有关，与扩散项的计算有关，
                                                                                                                 //在turbulenceProperties字典文件里的simulationType 设置为laminar时只有层流粘度
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);
//*************************************************************************************************// adj-fluid-a
Info<< "\nReading field pa\n" << endl;
volScalarField pa                                                                       // 定义伴随压力场pa
(
    IOobject
    (
        "pa",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
Info<< "\nReading field Ua\n" << endl;
volVectorField Ua                                                                       // 定义伴随速度场Ua
(
    IOobject
    (
        "Ua",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
Info<< "\nReading/calculating face flux field phi_a\n" << endl;
surfaceScalarField phia                                                              // 定义伴随通量场phi_a，面心标量场// 与adjointShapeOptimization一致
(
    IOobject
    (
        "phia",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(Ua) & mesh.Sf()                                       // 值为面上的速度u&面矢量sf
);
label paRefCell = 0;                                                                        // pa的参考压力
scalar paRefValue = 0.0;
setRefCell
(
    pa,
    simple.dict(),
    paRefCell,
    paRefValue
);
mesh.setFluxRequired(pa.name());

//*************************************************************************************************//
volScalarField fsensPowerDiss                                                                               // 定义木匾函数能量耗散对x的灵敏度，
(
    IOobject
    (
        "fsensPowerDiss",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -alphaMax*(1+qu)*qu/((qu+xh)*(qu+xh))*(U&Ua),                                    // 应为R^u(动量方程)中只有alpha_f*U项与xh有关，所以dJ/dxh=
                                                                                                                                             // ua&d(alpha_f*U)/dxh
    zeroGradientFvPatchScalarField::typeName
);
volScalarField fsenshPowerDiss(fsensPowerDiss);                                   // 定义能量耗散约束对xh的灵敏度，
