#include "AquaTopBGGenFilter.h"
#include <math.h>

/*==============================================================================
 * AquaTopBGGenFilter (Jun-Aug 2014 - Japan)
==============================================================================*/
AquaTopBGGenFilter::AquaTopBGGenFilter( TiXmlElement* _config, Filter* _input ): Filter( _config, _input )
{
	checkImage();
	invert = 0;
	adaptive = 0;
	resetOnInit = 1;
	if(useIntensityImage) background = new ShortImage( image->getWidth(), image->getHeight() );
	else background = new ShortImage( shortimage->getWidth(), shortimage->getHeight() );
	width = background->getWidth(); height = background->getHeight();
	maskimage = new ShortImage(width, height);
	tmpimage = new ShortImage(width, height);
	config->QueryIntAttribute( "Invert",   &invert   );
	config->QueryIntAttribute( "Adaptive", &adaptive );
	config->QueryIntAttribute( "MinSize",  &minsize  );
	config->QueryIntAttribute( "MaxSize",  &maxsize  );
	config->QueryIntAttribute( "PaperDepthDiff",  &paperdepthdiff);
	// setting variables for Configurator
	countOfOptions = 5; // quantity of variables that can be manipulated

	axis1 = (1,0,0);
	axis2 = (0,1,0);
	pos = (0,0,0);
}

AquaTopBGGenFilter::~AquaTopBGGenFilter() {
	delete background;
}

void AquaTopBGGenFilter::reset(int initialReset) {
	if( (initialReset == 1 && resetOnInit == 1) || initialReset == 0 ) {
		
		//TODO call function of other class here: 1) find rectangle (paper) (save coordinates of rectangle to send to client later), 2)extrapolate bg-values, 3)return new bg-image
		if(useIntensityImage) *background = *(input->getImage());
		else *background = *(input->getShortImage());
	}
}

int AquaTopBGGenFilter::process() {
	
	// only works for shortimage at this time (kinect uses shortimage)
	if(useIntensityImage)
	{
		image = input->getImage();
		return 0;
	}

	*tmpimage = *(input->getShortImage());
	unsigned short* tmpimagedata = tmpimage->getSData();

	// clone the input image to the IntensityImage
	tmpimage->convert(*image);

	//"threshold" data of background
	unsigned char* data = image->getData();
	for (int i = 0; i < width*height; i++) if(data[i] != 0) data[i] = 255;
	
	//local blob counter to differentiate blobs
	unsigned char value = 254;

	// Find the Paper blob(s) in the copy (use IntensityImage->integrate to get pos and size, overwrite values with a certain number) - see example in Bloblist.cc 138 and Blob.cc 24
	for (int i = 0; i < width*height; i++) if (data[i] == 255) try {
		// integrate the spot and fill it with the current counter value
		size = (int)image->integrate( Point(i%width,i/width), pos, axis1, axis2, 255, value );

		// if the spot is too small or too big, wipe it out and abort
		if ( (size < minsize) || ((maxsize != 0) && (size > maxsize)) ){
			image->integrate( Point(i%width,i/width), pos, axis1, axis2, value, 0 );
			continue;
		}
		// adjust counter
		value--;
		// did the frame-local blob counter overflow?
		if (value == 0) {
			value = 254;
			std::cerr << "Warning: too many blobs in AquaTopGenerator!" << std::endl;
		}
	} catch (...) { }
	
	//TODO temporary workaround. should be: if no paper blob found, just forward image (what if there are several papers? -.-)
	//if (value == 254)
		//return 0;

	//std::cout << "end-------------------------------------------------------------"<< std::endl;
	// use the imagecopy to find all pixels outside of the blobs in the original image and overwrite them with 0
	for (int i = 0; i < width*height; i++)	if(data[i] == 0) tmpimagedata[i] = 0;

	//copy shortimage to maskimage
	//*maskimage = *shortimage;

	paperblobcounter = 0;
	 //only keep the original values in the area of the blobs which are paper blobs, delete all paper blobs from mask image
	for(int j = 254; j > value; j--)
	{
		blobmin = 65535; blobmax = 0;
		for(int i = 0; i < width*height; i++) if(data[i] == j)
		{
			blobtmp = tmpimagedata[i];
			if(blobtmp < blobmin) blobmin = blobtmp;
			if(blobtmp > blobmax) blobmax = blobtmp;
		}
		//std::cout << "blobmin " << blobmin << " blobmax " << blobmax << std::endl;
		//std::cout << "blobdiff " << blobmax - blobmin << " for value " << static_cast<unsigned>(value) << std::endl;
		if(blobmax - blobmin < paperdepthdiff) //--> paper blob
		{
			paperblobcounter++;
		}
		else //--> non-paper blob
		{
			for (int i = 0; i < width*height; i++)	if(data[i] == j) tmpimagedata[i] = 0;
		}
	}

	
	// Generate BG for the whole surface (Least-Squares Plane fitting) if at least one paper blob was found
	if(paperblobcounter > 0)
		*background = *tmpimage;

	// BG subtraction
	*tmpimage = *(input->getShortImage());
	background->subtract( *tmpimage, *shortimage, invert );
	//if (adaptive) background->update( *(inputimg), *(maskimage) );
	return 0;
}

const char* AquaTopBGGenFilter::getOptionName(int option) {
	const char* OptionName = "";

	switch(option) {
	case 0:
		OptionName = "Invert";
		break;
	case 1:
		OptionName = "Adaptive";
		break;
	case 2:
		OptionName = "Minimum Size";
		break;
	case 3:
		OptionName = "Maximum Size";
		break;
	case 4:
		OptionName = "Paper Depth Diff";
		break;
	default:
		// leave OptionName empty
		break;
	}

	return OptionName;
}

double AquaTopBGGenFilter::getOptionValue(int option) {
	double OptionValue = -1.0;

	switch(option) {
	case 0:
		OptionValue = invert;
		break;
	case 1:
		OptionValue = adaptive;
		break;
	case 2:
		OptionValue = minsize;
		break;
	case 3:
		OptionValue = maxsize;
		break;
	case 4:
		OptionValue = paperdepthdiff;
		break;
	default:
		// leave OptionValue = -1.0
		break;
	}

	return OptionValue;
}

void AquaTopBGGenFilter::modifyOptionValue(double delta, bool overwrite) {
	switch(toggle) {
	case 0: // invert is a boolean value
		if(overwrite) {
			invert = (delta == 0 ? 0 : (delta == 1 ? 1 : invert));
		} else {
			invert += delta;
			invert = (invert < 0) ? 0 : (invert > 1) ? 1 : invert;
		}
		break;
	case 1: // adaptive is a boolean value
		if(overwrite) {
			adaptive = (delta == 0 ? 0 : (delta == 1 ? 1 : adaptive));
		} else {
			adaptive += delta;
			adaptive = (adaptive < 0) ? 0 : (adaptive > 1) ? 1 : adaptive;
		}
		break;
	case 2:
		if(overwrite) {
			minsize = (delta < 0) ? 0 : (delta > MAX_VALUE) ? MAX_VALUE : delta;
		} else {
			minsize += delta;
			minsize = (minsize < 0) ? 0 : (minsize > MAX_VALUE) ? MAX_VALUE : minsize;
		}
		break;
	case 3:
		if(overwrite) {
			maxsize = (delta < 0) ? 0 : (delta > MAX_VALUE) ? MAX_VALUE : delta;
		} else {
			maxsize += delta;
			maxsize = (maxsize < 0) ? 0 : (maxsize > MAX_VALUE) ? MAX_VALUE : maxsize;
		}
		break;
	case 4:
		if(overwrite) {
			paperdepthdiff = (delta < 0) ? 0 : (delta > MAX_VALUE) ? MAX_VALUE : delta;
		} else {
			paperdepthdiff += delta;
			paperdepthdiff = (paperdepthdiff < 0) ? 0 : (paperdepthdiff > MAX_VALUE) ? MAX_VALUE : paperdepthdiff;
		}
		break;
	}
}

TiXmlElement* AquaTopBGGenFilter::getXMLRepresentation() {
	
	TiXmlElement* XMLNode = new TiXmlElement( "AquaTopBGGenFilter" );
	
	XMLNode->SetAttribute( "Invert" , invert );
	XMLNode->SetAttribute( "Adaptive" , adaptive );
	XMLNode->SetAttribute( "MinSize",  minsize );
	XMLNode->SetAttribute( "MaxSize",  maxsize );
	XMLNode->SetAttribute( "PaperDepthDiff",  paperdepthdiff );
	
	return XMLNode;
}

// Fitting target: lowest sum of squared absolute error
// Fitting target value = 0.747158511837
//z = a + by + cx + dxy + f(x^2) + g(x^2)y
/*double AquaTopBGGenFilter::GetZValue(double x_in, double y_in)
{
	double temp;
	temp = 0.0;

	// coefficients
	double a = 1.8564803894488152E+03;
	double b = 5.2563093116817700E-02;
	double c = 1.2853574330037623E-01;
	double d = -2.7047119933664614E-04;
	double f = -1.9683753259302163E-04;
	double g = 4.7046643987525005E-07;

	temp += a;
	temp += b * y_in;
	temp += c * x_in;
	temp += d * x_in * y_in;
	temp += f * pow(x_in, 2.0);
	temp += g * pow(x_in, 2.0) * y_in;
	return temp;
}

//z = a + bx + cy
double AquaTopBGGenFilter::GetZValueLin(double x_in, double y_in)
{
	double temp;
	temp = 0.0;

	// coefficients
	double a = 1.8706125660880882E+03;
	double b = 1.9435859943767713E-02;
	double c = 1.6705938888233653E-02;

	temp = a + b * x_in + c * y_in;
	return temp;
}*/