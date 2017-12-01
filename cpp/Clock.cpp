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

  ==============================================================================
                      Implementation of class 'Clock'
  ==============================================================================*/

#include "Clock.h"
#include "GlobalsDefine.h"

// Constructor
Clock::Clock() 
{ 
	ct_ = 0;  
	st0_ = st1_ = false; 
}
  
// Start the clock//时钟工作并开启
void Clock::Start() 
{ 
	t0_ = clock();//返回从“开启这个程序进程”到“程序中调用clock()函数”时之间的CPU时钟计时单元（clock tick）数
	st0_ = true; 
}

// Stop the clock//时钟关闭
void Clock::Stop() 
{
	if (!st0_) //未调用
	{
		cerr << "\n*** Error *** In Clock :: Stop()";
		cerr << " : Method Start() must have been called before.\n";
	}

	if(!st1_) //已开启
	{
		t1_ = clock(); 
		ct_ += (double) (t1_ - t0_); //加上开启到现在的（clock tick）数
		st1_ = true; //关闭时钟
	}
}
  
// Resume the stopped clock //时钟开启
void Clock::Resume() 
{
	if (!st0_) 
	{
		cerr << "\n*** Error *** In Clock :: Resume()";
		cerr << " : Method Start() must have been called before.\n";
	}

	if (!st1_) {
        cerr << "\n*** Error *** In Clock::Resume()";
		cerr << " : Method Stop() must have been called before.\n";
	}
	else  
	{
		t0_ = clock();
		st1_ = false; //开启时钟
	}
}

// Clear the clock
void Clock::Clear() 
{ 
	ct_ = 0; 
	st0_ = st1_ = false;
}

// Return the elapsed time since the clock started 返回时钟工作后所开启的时间（秒数）
double Clock::ElapsedTime() 
{
	double elapsed;

	if (!st0_) { //时钟未工作
		cerr << "\n*** Error *** In Clock :: ElapsedTime()";
		cerr << " : Method Start() must have been called before.\n";
	}

	if (st1_)  // Timer has been stopped.
		elapsed = ct_;
	else
	{
		t1_ = clock(); 
		elapsed = ct_ + (double) (t1_ - t0_); 
	}
	return elapsed / CLOCKS_PER_SEC;// elapsed除以一秒钟内CPU运行的时钟周期数
}
