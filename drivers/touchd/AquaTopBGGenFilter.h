/*************************************************************************\
*    Part of the TISCH framework - see http://tisch.sourceforge.net/      *
*     Copyright (c) 2012 by Norbert Wiedermann, <wiederma@in.tum.de>      *
*   Licensed under GNU Lesser General Public License (LGPL) 3 or later    *
\*************************************************************************/

#ifndef _AQUATOPBGGENFILTER_H_
#define _AQUATOPBGGENFILTER_H_

#include <vector>
#include <iostream>

#include <GLUTWindow.h>
#include <TUIOOutStream.h>

#include "Filter.h"
#include "Blob.h"

using namespace std;

class AquaTopBGGenFilter: public Filter {
	public:
		AquaTopBGGenFilter( TiXmlElement* _config = 0, Filter* _input = 0 );
		virtual ~AquaTopBGGenFilter();
		virtual int process();
		virtual void reset(int initialReset);
		// Configurator
		virtual const char* getOptionName(int option);
		virtual double getOptionValue(int option);
		virtual void modifyOptionValue(double delta, bool overwrite);
		virtual TiXmlElement* getXMLRepresentation();
	protected:
		ShortImage* background;
		ShortImage* maskimage;
		ShortImage* tmpimage;
		int invert, adaptive;
		int width, height;

		::Vector axis1, axis2, pos;
		int size;						// blob size
		int minsize, maxsize;           // blob detection settings

		unsigned short blobmin, blobmax; //min and max depthvalue within a blob
		unsigned short blobtmp;

		int paperdepthdiff;
		int paperblobcounter;
};

#endif _AQUATOPBGGENFILTER_H_