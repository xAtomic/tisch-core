/*************************************************************************\
*    Part of the TISCH framework - see http://tisch.sourceforge.net/      *
*     Copyright (c) 2012 by Norbert Wiedermann, <wiederma@in.tum.de>      *
*   Licensed under GNU Lesser General Public License (LGPL) 3 or later    *
\*************************************************************************/

#ifndef _AQUATOPBGGENERATOR_H_
#define _AQUATOPBGGENERATOR_H_

#include "IntensityImage.h"
#include "ShortImage.h"
#include "RGBImage.h"


using namespace std;

class AquaTopBGGenerator{
public:
	AquaTopBGGenerator();
	virtual ~AquaTopBGGenerator();
	double GetZValue(double x_in, double y_in);
	double GetZValueLin(double x_in, double y_in);

protected:
	
	//IplImage* iplRGBImage;
	//IplImage* iplConverted;
    
};

#endif _AQUATOPBGGENERATOR_H_