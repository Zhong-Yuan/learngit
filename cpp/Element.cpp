/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.02, October 27, 2017                                        */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include <algorithm>
#include <Eigen/Eigen>
#include "Element.h"
#include "FourQ.h"
#include "GlobalsDefine.h"
#include <iostream>
#include <iomanip>

//! Virtual deconstructor
CElement::~CElement()
{
    if (!nodes)
        delete [] nodes;
    
    if (!ElementMaterial)
        delete [] ElementMaterial;

    if (!LocationMatrix)
        delete [] LocationMatrix;
}

//  Calculate the column height, used with the skyline storage scheme//计算单元对全局各自由度列高的贡献
//该单元位置矩阵各数和最小自由度数的差是该单元对全局各自由度的列高要求
void CElement::CalculateColumnHeight(unsigned int* ColumnHeight)
{	

//  Generate location matrix //计算单元的位置矩阵向量
    GenerateLocationMatrix();
//  Look for the row number of the first non-zero element
    unsigned int nfirstrow = INT_MAX; // <limits.h>中有INT_MAX，返回指定整数类型所能表示的最大值
    for (unsigned int i = 0; i < ND; i++)//返回位置矩阵向量中最小自由度号
        if (LocationMatrix[i] && LocationMatrix[i] < nfirstrow)
            nfirstrow = LocationMatrix[i];

//	Calculate the column height contributed by this element
	for (unsigned int i = 0; i < ND; i++)
	{
		unsigned int column = LocationMatrix[i];
		if (!column)
			continue;

		unsigned int Height = column - nfirstrow;
		if (ColumnHeight[column-1] < Height) ColumnHeight[column-1] = Height;
		//更新ColumnHeight
	}
	//cout << ColumnHeight[0] << setw(12) << ColumnHeight[1] <<setw(12) << ColumnHeight[2] <<endl;
}

#ifdef _EIGEN_ //使用Eigen
	//	Assemble the banded global stiffness matrix (skyline storage scheme)
	void CElement::assembly(CSkylineMatrix<double>* StiffnessMatrix)
	{
		//此函数不使用
	}
#else //不使用Eigen
	void CElement::assembly(double* Matrix, CSkylineMatrix<double>* StiffnessMatrix)
	{	
	//  Calculate element stiffness matrix //计算单元的刚度矩阵（不同种类单元函数不同）
		ElementStiffness(Matrix);
	//	Assemble global stiffness matrix
		for (unsigned int j = 0; j < ND; j++) 
		{
			unsigned int Lj = LocationMatrix[j];	// Global equation number corresponding to jth DOF of //the element单元j+1自由度对应全局号Lj
			if (!Lj) 
				continue;
	//  Address of diagonal element of column j in the one dimensional element stiffness matrix
		unsigned int DiagjElement = (j+1)*j/2 + 1;//单元的刚度矩阵的j对角元地址
			for (unsigned int i = 0; i <= j; i++)
			{
				unsigned int Li = LocationMatrix[i];	// Global equation number corresponding to ith DOF //of the element
				if (!Li) 
					continue;
				(*StiffnessMatrix)(Li,Lj) += Matrix[DiagjElement + j - i - 1];
			}
		}
		return;
	}
#endif
