/*************************************************************************\
*    Part of the TISCH framework - see http://tisch.sourceforge.net/      *
*  Copyright (c) 2006 - 2010 by Florian Echtler, TUM <echtler@in.tum.de>  *
*   Licensed under GNU Lesser General Public License (LGPL) 3 or later    *
\*************************************************************************/

#include "Filter.h"


BGSubFilter::BGSubFilter( TiXmlElement* _config, Filter* _input ): Filter( _config, _input ) {
	checkImage();
	background = new ShortImage( image->getWidth(), image->getHeight() );
	config->QueryIntAttribute( "Invert", &invert );
	config->QueryIntAttribute( "Adaptive", &adaptive );
}

BGSubFilter::~BGSubFilter() {
	delete background;
}

void BGSubFilter::link( Filter* _mask ) {
	mask = _mask;
}

void BGSubFilter::reset() {
	*background = *(input->getImage());
}

int BGSubFilter::process() {
	IntensityImage* inputimg = input->getImage();
	background->subtract( *(inputimg), *image, invert );
	if( adaptive ) background->update( *(inputimg), *(mask->getImage()) );
	result = background->intensity(); // FIXME: does 'invert' have to be factored in here?
	return 0;
}


// TODO: make hflip/vflip configurable
FlipFilter::FlipFilter( TiXmlElement* _config, Filter* _input ): Filter( _config, _input ) {
	checkImage();
}

// TODO: should be MMX-accelerated
int FlipFilter::process() {

	unsigned char* inbuf  = input->getImage()->getData();
	unsigned char* outbuf = image->getData();

	int width  = image->getWidth();
	int height = image->getHeight();

	int inoffset  = 0;
	int outoffset = width-1;

	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) outbuf[outoffset-x] = inbuf[inoffset+x];
		inoffset  += width;
		outoffset += width;
	}

	return 0;
}


// TODO: use result from bgsub filter for threshold adjustment
ThreshFilter::ThreshFilter( TiXmlElement* _config, Filter* _input ): Filter( _config, _input ) {
	checkImage();
	threshold_max = 128;
	threshold_min = 255;
	config->QueryIntAttribute( "Threshold_max", &threshold_max );
	config->QueryIntAttribute( "Threshold_min", &threshold_min );
}

int ThreshFilter::process() {
	input->getImage()->threshold( threshold_max, *image, threshold_min );
	return 0;
}


SpeckleFilter::SpeckleFilter( TiXmlElement* _config, Filter* _input ): Filter( _config, _input ) {
	checkImage();
	noiselevel = 7;
	config->QueryIntAttribute( "NoiseLevel", &noiselevel );
}

int SpeckleFilter::process() {
	input->getImage()->despeckle( *image, noiselevel );
	return 0;
}


SplitFilter::SplitFilter( TiXmlElement* _config, Filter* _input ): Filter( _config, _input ) {
	checkImage();
	image2 = NULL;
	reset();
}

void SplitFilter::reset() {
	incount = outcount = 0;
}	

int SplitFilter::process() {
	incount++;
	if (incount % 2) {
		*image = *(input->getImage());
		return 1;
	} else {
		image2 = input->getImage();
		return 0;
	}
}

// TODO: add intensity heuristic (e.g. 2nd image is the brighter one)
IntensityImage* SplitFilter::getImage() {
	outcount++;
	if (outcount % 2) return image;
	else return image2 ? image2 : image;
}
