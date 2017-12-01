/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.02, October 27, 2017                                        */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#pragma once

#include "Outputter.h"

#include <string>
#include <fstream>

using namespace std;

//! Class LoadData is used to store load data//载荷类
class CLoadCaseData
{
public:

	unsigned int nloads;	//!< Number of concentrated loads in this load case 集中载荷个数
	unsigned int* node;		//!< Node number to which this load is applied 作用节点
	unsigned int* dof;		//!< Degree of freedom number for this load component 载荷分量对应的自由度
	double* load;			//!< Magnitude of load 载荷大小
	double grav;//重力加速度
	int Dgrav;//重力加速度所在的方向 1:x y z -x -y 6:-z

public:

	CLoadCaseData() : nloads(0), node(NULL), dof(NULL), load(NULL) ,grav(0),Dgrav(1){};//默认重力朝x方向
	~CLoadCaseData();

//!	Set nloads, and new array node, dof and load
	void Allocate(unsigned int num);

//!	Read load case data from stream Input
	bool Read(ifstream& Input, unsigned int lcase);

//!	Write load case data to stream
	void Write(COutputter& output, unsigned int lcase);
};
