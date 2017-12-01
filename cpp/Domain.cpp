/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.02, October 27, 2017                                        */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Domain.h"
#include "Bar.h"
#include "FourQ.h"
#include "Material.h"
#include "GlobalsDefine.h"

#include <iomanip>
#include <iostream>
#include <Eigen/Eigenvalues>

using namespace std;

//	Clear an array
template <class type> void clear( type* a, unsigned int N )
{
	for (unsigned int i = 0; i < N; i++)
		a[i] = 0;
}

CDomain* CDomain::_instance = nullptr;
//使用NULL实际上是定义为0,使用nullptr更好,能避免 0 指针的二义性的问题

//	Constructor
CDomain::CDomain()
{
	Title[0] = '0';
	MODEX = 0;

	NUMNP = 0;
	NodeList = nullptr;
	
	NUMEG = 0;
	EleGrpList = nullptr;

	NLCASE = 0;
	NLOAD = nullptr;
	LoadCases = nullptr;
	
	NEQ = 0;
	NWK = 0;
	MK = 0;

	Force = nullptr; 
	StiffnessMatrix = nullptr;
//	ElementTypes = nullptr;
//	NUME = NULL;
//	ElementSetList = NULL;
//	
//	NUMMAT = NULL;
//	MaterialSetList = NULL;
//	
}

//	Desconstructor
CDomain::~CDomain()
{
	delete [] NodeList;

	delete [] EleGrpList;

	delete [] NLOAD;
	delete [] LoadCases;

	delete [] Force;
	delete StiffnessMatrix;
//	delete [] ElementTypes;
//	delete [] NUME;
//	delete [] ElementSetList;
//
//	delete [] NUMMAT;
//	delete [] MaterialSetList;
//
}

//	Return pointer to the instance of the Domain class
CDomain* CDomain::Instance()
{
	if (!_instance) 
		_instance = new CDomain();
	
	return _instance;
}

//	Read domain data from the input data file
bool CDomain::ReadData(string FileName, string OutFile)
{
	Input.open(FileName);
	if (!Input) 
	{
		cerr << "*** Error *** File " << FileName << " does not exist !" << endl;
		exit(3);
	}

	COutputter* Output = COutputter::Instance(OutFile);

//	Read the heading line
	Input.getline(Title, 256);//用于读取一行字符直到换行符
	Output->OutputHeading();//输出Title和时间

//	Read the control line //读取节点总数 单元组数量 载荷数 求解模式
	Input >> NUMNP >> NUMEG >> NLCASE >> MODEX;

//	Read nodal point data
	if (ReadNodalPoints())//建立节点列表 NodeList = new CNode[NUMNP];
        Output->OutputNodeInfo();//输出节点信息
    else
        return false;

//	Update equation number
	CalculateEquationNumber();//将NodeList中的bcode改为全局方程号
	Output->OutputEquationNumber();

//	Read load data
	if (ReadLoadCases())//建立载荷类列表LoadCases并赋值 LoadCases = new CLoadCaseData[NLCASE];
        Output->OutputLoadInfo();
    else
        return false;

//	Read element data 读取单元组数据 EleGrpList = new CElementGroup[NUMEG];
	if (ReadElements())
        Output->OutputElementInfo();
    else
        return false;
	return true;
}

//	Read nodal point data //建立NodeList
bool CDomain::ReadNodalPoints()
{

//	Read nodal point data lines
	NodeList = new CNode[NUMNP];

//	Loop over for all nodal points
	for (unsigned int np = 0; np < NUMNP; np++)
		if (!NodeList[np].Read(Input, np))
			return false;

	return true;
}

//	Calculate global equation numbers corresponding to every degree of freedom of each node
void CDomain::CalculateEquationNumber()
{
	NEQ = 0;
	for (unsigned int np = 0; np < NUMNP; np++)	// Loop over for all node
	{
		for (unsigned int dof = 0; dof < CNode::NDF; dof++)	// Loop over for DOFs of node np
		{
			if (NodeList[np].bcode[dof]) //0：自由 1：固定
				NodeList[np].bcode[dof] = 0;
			else
			{
				NEQ++;
				NodeList[np].bcode[dof] = NEQ;
			}
		}
	}
}

//	Read load case data
bool CDomain::ReadLoadCases()
{
//	Read load data lines
	LoadCases = new CLoadCaseData[NLCASE];	// List all load cases 建立载荷类列表LoadCases

//	Loop over for all load cases
	for (unsigned int lcase = 0; lcase < NLCASE; lcase++) 
		if (!LoadCases[lcase].Read(Input, lcase))
			return false;

	return true;
}

// Read element data
bool CDomain::ReadElements()
{
    EleGrpList = new CElementGroup[NUMEG];

//	Loop over for all element group
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
        if (!EleGrpList[EleGrp].Read(Input))
            return false;
    
    return true;
}
////读取单元数据放于MaterialSetList材料设置列表中
//bool CDomain::ReadElements()
//{
////	Read element group control line NUMEG 单元组数量
//	ElementTypes = new unsigned int[NUMEG];	// Element type of each group 单元组对应的单元种类
//	NUME = new unsigned int[NUMEG];			// Number of elements in each group 单元组的单元数
//	ElementSetList = new CElement*[NUMEG];	// Element list in each group 单元设置列表 每行表示一个单元组设置
//
//	NUMMAT = new unsigned int[NUMEG];		// Material set number of each group 单元组对应的材料集编号
//	MaterialSetList = new CMaterial*[NUMEG];// Material list in each group 材料设置列表 每行表示一个单元组材料设置
//
////	Loop over for all element group 不断读取单元组
//	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
//	{
//		Input >> ElementTypes[EleGrp] >> NUME[EleGrp] >> NUMMAT[EleGrp];
//		//单元组对应的单元种类 单元组的单元数 单元组对应的材料集编号
//		switch (ElementTypes[EleGrp])
//		{
//		case 1:	// Bar element
//			if (!ReadBarElementData(EleGrp))//读取第EleGrp单元组的杆单元数据：杆材料和单元
//                return false;
//			break;
//		default:	// Invalid element type
//            cout << "*** Error *** Elment type " << ElementTypes[EleGrp] << " of group "
//                 << EleGrp+1 << " has not been implemented.\n\n";
//			return false;
//		}
//	}
//	return true;
//}
//
////	Read bar element data from the input data file //读取第EleGrp单元组的杆单元数据：杆材料和单元
//bool CDomain::ReadBarElementData(unsigned int EleGrp)
//{
////	Read material/section property lines
//	MaterialSetList[EleGrp] = new CBarMaterial[NUMMAT[EleGrp]];	// Materials for group EleGrp
//	//NUMMAT是单元组中的材料数目
//	//设置第EleGrp单元组的材料列表
//    CBarMaterial* mlist = (CBarMaterial*) MaterialSetList[EleGrp];
//	//取别名
////	Loop over for all material property sets in group EleGrp
//	for (unsigned int mset = 0; mset < NUMMAT[EleGrp]; mset++)
//		if (!mlist[mset].Read(Input, mset))//读取杆材料的杨氏模量和截面积
//			return false;
//
////	Read element data lines
//	ElementSetList[EleGrp] = new CBar[NUME[EleGrp]];// Elements of group EleGrp
//	//创建杆单元，定义单元列表的EleGrp行
////	Loop over for all elements in group EleGrp
//	for (unsigned int Ele = 0; Ele < NUME[EleGrp]; Ele++)
//		if (!ElementSetList[EleGrp][Ele].Read(Input, Ele, MaterialSetList[EleGrp], NodeList))//读取杆单元
//			return false;
//	return true;
//}

//	Calculate column heights 计算各自由度的列高ColumnHeights和最大列高（带宽）（包括对角元）
void CDomain::CalculateColumnHeights()
{
	unsigned int* ColumnHeights = StiffnessMatrix->GetColumnHeights();//返回列高数组
//clear(ColumnHeights, NEQ);	// Set all elements to zero 清空向量

	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)		//	Loop over for all element groups
	{
        CElementGroup* ElemengGrp = &EleGrpList[EleGrp];
        unsigned int NUME = ElemengGrp->GetNUME();
		//这里区分单元的种类!!!!!!!!!!!!!!!!
		unsigned int ElementType = ElemengGrp->GetElementType();
		switch (ElementType) //这种方式是否可以优化??
		{
			case 1:
				CalculateColumnHeightForCBar();
				break;
			case 2:
				CalculateColumnHeightForFourQ();
				break;
		}
		
	}
//	Maximum half bandwidth ( = max(ColumnHeights) + 1 )
//计算最大列高（包括对角元）
	MK = ColumnHeights[0];

	for (unsigned int i=1; i<NEQ; i++)
		if (MK < ColumnHeights[i])
			MK = ColumnHeights[i];
	MK = MK + 1;
#ifdef _DEBUG_
	COutputter* Output = COutputter::Instance();//COutputter类已经在ReadData时建立
	Output->PrintColumnHeights();
#endif
}

void CDomain::CalculateColumnHeightForCBar()
{
	CBar* ElementList = dynamic_cast<CBar*> (ElemengGrp->GetElementList());
	for (unsigned int Ele = 0; Ele < NUME; Ele++)	//	Loop over for all elements in group EleGrp
	ElementList[Ele].CalculateColumnHeight(ColumnHeights); //计算单元组中的一个单元对列高的贡献//该单元位置矩阵各数和最小自由度数的差是该单元对全局各自由度的列高要求
}

void CDomain::CalculateColumnHeightForFourQ()
{
	FourQ* ElementList = dynamic_cast<FourQ*> (ElemengGrp->GetElementList());
	for (unsigned int Ele = 0; Ele < NUME; Ele++)	//	Loop over for all elements igroup EleGrp
	ElementList[Ele].CalculateColumnHeight(ColumnHeights); //计算单元组中的一个单元高的贡献//该单元位置矩阵各数和最小自由度数的差是该单元对全局各自由度的列高要求
	break;
}

//	Calculate address of diagonal elements in banded matrix
//	Caution: Address is numbered from 1 !
void CDomain::CalculateDiagnoalAddress()
{
    unsigned int* ColumnHeights = StiffnessMatrix->GetColumnHeights();
    unsigned int* DiagonalAddress = StiffnessMatrix->GetDiagonalAddress();
//	clear(DiagonalAddress, NEQ + 1);	// Set all elements to zero

//	Calculate the address of diagonal elements
//	M(0) = 1;  M(i+1) = M(i) + H(i) + 1 (i = 0:NEQ)
	DiagonalAddress[0] = 1;
	for (unsigned int col = 1; col <= NEQ; col++)
		DiagonalAddress[col] = DiagonalAddress[col - 1] + ColumnHeights[col-1] + 1;

//	Number of elements in banded global stiffness matrix
	NWK = DiagonalAddress[NEQ] - DiagonalAddress[0];

#ifdef _DEBUG_
	COutputter* Output = COutputter::Instance();
	Output->PrintDiagonalAddress();
#endif

}

//	Assemble the banded gloabl stiffness matrix
void CDomain::AssembleStiffnessMatrix()
{
//	Loop over for all element groups
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{	
		//unsigned int size = ElementSetList[EleGrp][0].SizeOfStiffnessMatrix();
		//对3维杆单元的刚度阵有21自由度
		CElementGroup* ElemengGrp = &EleGrpList[EleGrp];
        unsigned int NUME = ElemengGrp->GetNUME();
		unsigned int ElementType = ElemengGrp->GetElementType();
		switch (ElementType) //这里是否可以
		{
			case 1:
				AssembleStiffnessMatrixForCBar();
				break;
			case 2:
				AssembleStiffnessMatrixForFourQ();
				break;
		}
	}
#ifdef _DEBUG_
	COutputter* Output = COutputter::Instance();
	Output->PrintStiffnessMatrix();
#endif
}

void CDomain::AssembleStiffnessMatrixForCBar()
{
	CBar* ElementList = dynamic_cast<CBar*> (ElemengGrp->GetElementList());
	#ifdef _EIGEN_ //使用Eigen
		for (unsigned int Ele = 0; Ele < NUME; Ele++)
			//将各单元组装至全局刚度阵
			ElementList[Ele].ElementStiffness(StiffnessMatrix);//这里直接进入CBar的ElementStiffness函数，将组装工作也分散至各单元中! 可优化
	#else //不使用Eigen
		//Loop over for all elements in group EleGrp
		unsigned int size = ElementList[0].SizeOfStiffnessMatrix();
		double* Matrix = new double[size]; //尽量减少堆申请的次数
		for (unsigned int Ele = 0; Ele < NUME; Ele++)
			//将各单元组装至全局刚度阵
			ElementList[Ele].assembly(Matrix, StiffnessMatrix);
		delete [] Matrix;
	#endif
}

void CDomain::AssembleStiffnessMatrixForFourQ()
{
	FourQ* ElementList = dynamic_cast<FourQ*> (ElemengGrp->GetElementList());
	#ifdef _EIGEN_ //使用Eigen
		for (unsigned int Ele = 0; Ele < NUME; Ele++)
			//将各单元组装至全局刚度阵
			ElementList[Ele].assembly(StiffnessMatrix);
		//Eigen::Matrix<double, size, size>* Matrix = new Eigen::Matrix<double, size, size>();//临时存储矩阵
	#else //不使用Eigen
		unsigned int size = ElementList[0].SizeOfStiffnessMatrix();
		double* Matrix = new double[size];
		//Loop over for all elements in group EleGrp
		for (unsigned int Ele = 0; Ele < NUME; Ele++)
			//将各单元组装至全局刚度阵
			ElementList[Ele].assembly(Matrix, StiffnessMatrix);
		delete [] Matrix;
	#endif
}

bool CDomain::BarDisGForce(unsigned int EleGrp,double g,int Dg)
{	
	//g是重力加速度 Dg是作用方向 1:6
	cout<<EleGrp<<g<<Dg<<endl;//
	//CBarMaterial* mlist =(CBarMaterial*) MaterialSetList[EleGrp];//取别名
	CBar* blist = dynamic_cast<CBar *>(EleGrpList[EleGrp].GetElementList());
		
	double gravForce = 0;//每根杆的重力
	int sign;
	int D;//xyz 012
	unsigned int dof;
	D = (Dg-1) % 3;
	sign = 1;//重力作用的正方向
	if ( Dg>=4 ){ sign=-1; }
	
	int NUME = EleGrpList[EleGrp].GetNUME();
	for (unsigned int Ele = 0; Ele < NUME; Ele++)
	{	
		CBarMaterial* CBarM =dynamic_cast<CBarMaterial *>(blist[Ele].ElementMaterial);
		gravForce = CBarM->Rho*g*blist[Ele].GetLength()*CBarM->Area;
		//将该杆重力分散到节点上
		gravForce = gravForce/2;
		//左节点对应的方程号
		dof = blist[Ele].nodes[0]->bcode[D]; //注意设置方向
		if(dof !=0) { Force[dof - 1] += sign*gravForce; }
		//右节点对应的方程号
		dof = blist[Ele].nodes[1]->bcode[D];
		if(dof !=0) { Force[dof - 1] += sign*gravForce; }
	}	
	return true;
}

bool CDomain::DisGForceall(double g,int Dg)//为所有单元进行重力分配
{	
	for(unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)//对单元组循环
	{	
		switch (EleGrpList[EleGrp].GetElementType())
		{
		case 1:	// Bar element
			if (!BarDisGForce(EleGrp,g,Dg))//对Bar的重力进行分配
                return false;
			break;
		}
	}
	return true;
}
//	Assemble the global nodal force vector for load case LoadCase //对第LoadCase载荷类组装力矢量
bool CDomain::AssembleForce(unsigned int LoadCase)
{
	if (LoadCase > NLCASE) 
		return false;
	CLoadCaseData* LoadData = &LoadCases[LoadCase - 1];
	clear(Force, NEQ);//清空力向量
	//Loop over for all concentrated loads in load case LoadCase
	//对集中载荷个数循环
	for (unsigned int lnum = 0; lnum < LoadData->nloads; lnum++)
	{
		//载荷对应的方程号
		unsigned int dof = NodeList[LoadData->node[lnum] - 1].bcode[LoadData->dof[lnum] - 1];
        if(dof) // The DOF is activated
			Force[dof - 1] += LoadData->load[lnum];
	}
	//将重力分散至各个节点
	if(LoadCases[LoadCase - 1].grav != 0)//存在重力加速度
	{
		if(!DisGForceall(LoadCases[LoadCase - 1].grav,LoadCases[LoadCase - 1].Dgrav))
			return false;
	}
	
	return true;
}

//	Allocate storage for matrices Force, ColumnHeights, DiagonalAddress and StiffnessMatrix
//	and calculate the column heights and address of diagonal elements
//为矩阵力，列高，对角元地址和刚度矩阵分配空间 ???修正??????
//计算列高度和对角线单元的地址 
void CDomain::AllocateMatrices()
{
//	Allocate for global force/displacement vector 全局力/位移矢量
	Force = new double[NEQ];
	clear(Force, NEQ);

//  Create the banded stiffness matrix //初始化列高和对角地址矩阵
    StiffnessMatrix = new CSkylineMatrix<double>(NEQ);

//	Calculate column heights//计算列高
	CalculateColumnHeights();

//	Calculate address of diagonal elements in banded matrix
	CalculateDiagnoalAddress();

//	Allocate for banded global stiffness matrix //计算出列高和分配矩阵存储空间
	StiffnessMatrix->Allocate();

	COutputter* Output = COutputter::Instance();
	Output->OutputTotalSystemData();
}
