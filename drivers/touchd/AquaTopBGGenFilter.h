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
		void send( TUIOOutStream* tuio );
		virtual void reset(int initialReset);
		// Configurator
		virtual const char* getOptionName(int option);
		virtual double getOptionValue(int option);
		virtual void modifyOptionValue(double delta, bool overwrite);
		virtual TiXmlElement* getXMLRepresentation();

	protected:

		struct Matrix
		{
			long double m00;
			long double m10;
			long double m20;
			long double m01;
			long double m11;
			long double m21;
			long double m02;
			long double m12;
			long double m22;
		};

		void GenerateBackground();

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
		bool hasnonpaperblob;

		std::vector<Blob>* paperblobs;
		std::vector<Blob>* tmpblobs;
		std::vector<::Vector>* papercorners;

		int id;

		IntensityImage* testimage1;
		IntensityImage* testimage2;
		IntensityImage* testimage3;
		IntensityImage* testimage4;
		IntensityImage* testimage5;

		int houghTr;
		int houghTtheta;
		int houghTanglediff;
};

#endif _AQUATOPBGGENFILTER_H_