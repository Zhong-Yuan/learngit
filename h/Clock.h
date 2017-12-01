/*==============================================================================
                                    MPM3D++
           C++ code for Three Dimensional Material Point Method 
  ==============================================================================

   Copyright (C) 2006 - 

   Computational Dynamics Group
   Department of Engineering Mechanics
   Tsinghua University
   Beijing 100084, P. R. China

   Email: xzhang@tsinghua.edu.cn
*/

#pragma once

#include <time.h>
#include <iostream>

using namespace std;  

//! Clock class for timing 时钟：工作：开启/停止
class Clock
{

private:

	clock_t t0_, t1_; //开启时的时刻 停止时的时刻
	double ct_; //返回开启时间的累计
	bool st0_;   //!< Flag for Start method 时钟是否工作
	bool st1_;   //!< Flag for Stop method 时钟是否开启

public:

//!	Constructor
	Clock();
  
//!	Start the clock
	void Start();

//!	Stop the clock
	void Stop();
  
//!	Resume the stoped clock
	void Resume();

//!	Clear the clock
	void Clear();

//!	Return the elapsed time since the clock started
	double ElapsedTime();

};
