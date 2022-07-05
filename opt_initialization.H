/********************************************************************************/
/*****************必要变量声明及初始化,如MMA***********************************/
/********************************************************************************/


Info<< "\nClear output files of last solve \n" << endl;
std::remove("Voluse.txt");
std::remove("Time.txt");
std::remove("PowerDissRatio.txt"); 
std::remove("Change.txt");

int i = 0, opt = 0;                     // opt用于优化过程参数的放松，在"update.H"里opt++,这里先进行初始化
double del = 0;                       // del就是heviside过滤里的惩罚强度参数betta
scalar n_dof = xh.size();               //变量个数，不考虑并行
scalar m_con = 1;                           // 约束个数
scalar eta0, eta5(0.5), eta1, y0, y5;//eta5就是heviside过滤的截断值eta_h，小于该值的伪密度会倾向于罚为0；大于该值的倾向于罚为1

// n = xh.size(); // n:length of x in the i-th core
// N = xh.size(); // N:design freedom;
// reduce(N, sumOp<scalar>());// 把所有核中的值累加在一起，reduce()函数中sumOp描述累加，还有minOp，maxOp等

double PowerDiss(0), V(0), area = gSum(mesh.V());                            //定义了柔度C，体积V，area是单元体积和，gSum() sums over all the processors in a parallel run
double len = std::pow(area / n_dof, 1.0 / 3);                            // 这个len计算出所有单元的平均尺寸
Info << "\n网格的平均尺寸为：" << len << "\n"<<endl;
b.value() = 1.0 / (coff_PDE * len) / (coff_PDE * len); // b就是PDE求导数中的1/( r^2 ),coff_PDE * len=r,r为coff_PDE倍的单元长度;

/*****************************MMA相关*******************************************/

Info<< "Initialization of MMA \n" << endl;

// Math helpers
double Min(double d1, double d2);
double Max(double d1, double d2);
int Min(int d1, int d2);
int Max(int d1, int d2);
double Abs(double d1);
scalar diff(volScalarField gamma, scalarField V, double del, double eta, int n);

int n = n_dof;
int m = m_con;

double Xmin = 0.0;
double Xmax = 1.0;
double movlim = 0.2;

double *xmma = new double[n];
double *xold = new double[n];
double *df = new double[n];
double *g = new double[m];
double *dg = new double[m * n];
double *xmin = new double[n];
double *xmax = new double[n];

for(i=0;i<n;i++)
{
    xmma[i]=x[i];
}
  	// Initialize MMA
	MMASolver *mma = new MMASolver(n,m);
//------------------------------------------------------------------------------------------------------------//


