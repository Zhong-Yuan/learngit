/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.02, October 27, 2017                                        */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Material.h"
#include "GlobalsDefine.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <Eigen/Eigen>


using namespace std;

//	Read material data from stream Input//读取杆材料数据
bool CBarMaterial::Read(ifstream& Input, unsigned int mset)
{
	Input >> nset;	// Number of property set 杆材料编号

	if (nset != mset + 1)
	{
		cerr << "*** Error *** Material sets must be inputted in order !" << endl 
			 << "    Expected set : " << mset + 1 << endl
			 << "    Provided set : " << nset << endl;

		return false;
	}
	Input >> E >> Area >> Rho;	// Young's modulus and section area 材料密度 
	return true;
}

//	Write material data to Stream
void CBarMaterial::Write(COutputter& output, unsigned int mset)
{
	output << setw(5) << mset+1 << setw(16) << E << setw(16) << Area << setw(16) << Rho <<endl;
}

//	Read material data from stream Input//读取4Q材料数据
bool FourQMaterial::Read(ifstream& Input, unsigned int mset)
{
	Input >> nset;	// Number of property set 4Q材料编号

	if (nset != mset + 1)
	{
		cerr << "*** Error *** Material sets must be inputted in order !" << endl 
			 << "    Expected set : " << mset + 1 << endl
			 << "    Provided set : " << nset << endl;
		return false;
	}
	Input >> E >> nu >> Rho;	// 杨氏模量 泊松比 材料密度
	CalcD();//计算D矩阵
	return true;
}

//	Write material data to Stream
void FourQMaterial::Write(COutputter& output, unsigned int mset)
{
	output << setw(5) << mset+1 << setw(16) << E << setw(16) << nu << setw(16) << Rho <<endl;
}

void FourQMaterial::CalcD()
{
	D <<  1, nu, 0,
		 nu,  1, 0,
		  0,  0, (1-nu)/2;
	D*=E/(1-nu*nu);
}





