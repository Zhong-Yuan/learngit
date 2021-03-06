/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.02, October 27, 2017                                        */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Bar.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <Eigen/Eigen>
#include "GlobalsDefine.h"

using namespace std;

//	Constructor
CBar::CBar()
{
	NEN = 2;	// Each element has 2 nodes
	nodes = new CNode*[NEN];
    
    ND = 6;
    LocationMatrix = new unsigned int[ND];

	ElementMaterial = nullptr;
	Len = 0;    // 0长度是非法值，调用时可以检测
}

//	Desconstructor
CBar::~CBar()
{
	//不需要显式地调用基类的析构函数，系统会自动隐式调用。
	//先执行派生类析构函数的函数体，再调用基类的析构函数。
}

//	Read element data from stream Input 读取并配置单元组中第Ele号杆单元的材料配置号（在对应材料组中的编号），左右节点号
bool CBar::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList)
//MaterialSets实参是CBarMaterial* MaterialSetList[EleGrp]
{
	unsigned int N;
	Input >> N;	// element number
	if (N != Ele + 1)
	{
		cerr << "*** Error *** Elements must be inputted in order !" << endl 
			 << "    Expected element : " << Ele + 1 << endl
			 << "    Provided element : " << N << endl;
		//cerr：输出到标准错误的ostream对象，常用于程序错误信息
		return false;
	}
	unsigned int MSet;	// Material property set number
	unsigned int N1, N2;	// Left node number and right node number

	Input >> N1 >> N2 >> MSet;
	//ElementMaterial = &(dynamic_cast<CBarMaterial*>(MaterialSets))[MSet - 1];
	//CBarMaterial* ElementMaterial //dynamic_cast 将一个基类对象指针（或引用）cast到继承类指针
	ElementMaterial = dynamic_cast<CBarMaterial*>(MaterialSets) + MSet - 1; //对指针加法
	nodes[0] = &NodeList[N1 - 1];
	nodes[1] = &NodeList[N2 - 1];

	return true;
}

//	Write element data to stream
void CBar::Write(COutputter& output, unsigned int Ele)
{
	output << setw(5) << Ele+1 << setw(11) << nodes[0]->NodeNumber 
		   << setw(9) << nodes[1]->NodeNumber << setw(12) << ElementMaterial->nset << endl;
}


//!获得杆的长度
double CBar::GetLength()
{	
	return Len;
}
//计算杆的长度 //在CBar::ElementStiffness中已经对长度进行计算!
double CBar::CalcLength()
{
	double Dx;
	double Dy;
	double Dz;
	double Len;
	Dx = nodes[0]->XYZ[0]-nodes[1]->XYZ[0];
	Dy = nodes[0]->XYZ[1]-nodes[1]->XYZ[1];
	Dz = nodes[0]->XYZ[2]-nodes[1]->XYZ[2];
	Len = sqrt(Dx*Dx + Dy*Dy + Dz*Dz);	
}

//  Generate location matrix: the global equation number that corresponding to each DOF of the element
//	Caution:  Equation number is numbered from 1 !
void CBar::GenerateLocationMatrix()
{
    unsigned int i = 0;
    for (unsigned int N = 0; N < NEN; N++) //每个单元节点数循环
        for (unsigned int D = 0; D < 3; D++) //每个节点自由度循环
            LocationMatrix[i++] = nodes[N]->bcode[D];//这时bcode存储的是全局自由度编号
}


//	Return the size of the element stiffness matrix (stored as an array column by column)
//	For 2 node bar element, element stiffness is a 6x6 matrix, whose upper triangular part
//	has 21 elements
unsigned int CBar::SizeOfStiffnessMatrix() 
	{ 
		#ifdef _EIGEN_ //使用Eigen
			return 6;
		#else //不使用Eigen
			return 21; 
		#endif
	}


#ifdef _EIGEN_ //使用Eigen
//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CBar::ElementStiffness(CSkylineMatrix<double>* StiffnessMatrix)
{
	//待实现
}
#else //不使用Eigen
//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CBar::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());

//	Calculate bar length
	double DX[3];		//	dx = x2-x1, dy = y2-y1, dz = z2-z1
	for (unsigned int i = 0; i < 3; i++)
		DX[i] = nodes[1]->XYZ[i] - nodes[0]->XYZ[i];

	double DX2[6];	//  Quadratic polynomial (dx^2, dy^2, dz^2, dx*dy, dy*dz, dx*dz)
	DX2[0] = DX[0] * DX[0];
	DX2[1] = DX[1] * DX[1];
	DX2[2] = DX[2] * DX[2];
	DX2[3] = DX[0] * DX[1];
	DX2[4] = DX[1] * DX[2];
	DX2[5] = DX[0] * DX[2];

	double L2 = DX2[0] + DX2[1] + DX2[2];
	Len = sqrt(L2);
//	Calculate element stiffness matrix
	CBarMaterial* material = dynamic_cast<CBarMaterial*>(ElementMaterial);	// Pointer to material of the element
	double k = material->E * material->Area / Len / L2;
	Matrix[0] = k*DX2[0];
	Matrix[1] = k*DX2[1];
	Matrix[2] = k*DX2[3];
	Matrix[3] = k*DX2[2];
	Matrix[4] = k*DX2[4];
	Matrix[5] = k*DX2[5];
	Matrix[6] = k*DX2[0];
	Matrix[7] = -k*DX2[5];
	Matrix[8] = -k*DX2[3];
	Matrix[9] = -k*DX2[0];
	Matrix[10] = k*DX2[1];
	Matrix[11] = k*DX2[3];
	Matrix[12] = -k*DX2[4];
	Matrix[13] = -k*DX2[1];
	Matrix[14] = -k*DX2[3];
	Matrix[15] = k*DX2[2];
	Matrix[16] = k*DX2[4];
	Matrix[17] = k*DX2[5];
	Matrix[18] = -k*DX2[2];
	Matrix[19] = -k*DX2[4];
	Matrix[20] = -k*DX2[5];
}
#endif

//	Calculate element stress //杆单元的应力计算
void CBar::ElementStress(double* stress, double* Displacement)
{
	CBarMaterial* material = dynamic_cast<CBarMaterial*>(ElementMaterial);	// Pointer to material of the element

	double DX[3];	//	dx = x2-x1, dy = y2-y1, dz = z2-z1
	double L2 = 0;	//	Square of bar length (L^2)

	for (unsigned int i = 0; i < 3; i++)
	{
		DX[i] = nodes[1]->XYZ[i] - nodes[0]->XYZ[i];
		L2 = L2 + DX[i]*DX[i];
	}

	double S[6];
	for (unsigned int i = 0; i < 3; i++)
	{
		S[i] = -DX[i] * material->E / L2;
		S[i+3] = -S[i];
	}
	
	*stress = 0.0;
	for (unsigned int i = 0; i < 6; i++)
	{
		if (LocationMatrix[i])
			*stress += S[i] * Displacement[LocationMatrix[i]-1];
	}
}
