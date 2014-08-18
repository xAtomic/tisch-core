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
	config->QueryIntAttribute( "HoughTr",  &houghTr);
	config->QueryIntAttribute( "HoughTtheta",  &houghTtheta);
	config->QueryIntAttribute( "HoughTanglediff",  &houghTanglediff);
	// setting variables for Configurator
	countOfOptions = 8; // quantity of variables that can be manipulated

	axis1 = (1,0,0);
	axis2 = (0,1,0);
	pos = (0,0,0);
	id = 0;
	hasnonpaperblob = false;

	tmpblobs = new std::vector<Blob>;
	paperblobs = new std::vector<Blob>;
	papercorners = new std::vector<::Vector>;

				testimage1 = new IntensityImage(*input->getImage());
			testimage2 = new IntensityImage(*input->getImage());
			testimage3 = new IntensityImage(0.15*testimage2->getWidth(), 0.5*testimage2->getHeight());
			testimage4 = new IntensityImage(0.15*testimage2->getWidth(), 0.5*testimage2->getHeight());
			testimage5 = new IntensityImage(0.15*testimage2->getWidth(), 0.5*testimage2->getHeight());
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
	rgbimage = input->getRGBImage(); // get pointer from previous filter, do nothing
	
	// only works for shortimage at this time (kinect uses shortimage)
	if(useIntensityImage)
	{
		image = input->getImage();
		return 0;
	}
	

	// clear lists
	tmpblobs->clear();
	paperblobs->clear();
	papercorners->clear();

	*tmpimage = *(input->getShortImage());
	unsigned short* tmpimagedata = tmpimage->getSData();

	// clone the input image to the IntensityImage
	tmpimage->convert(*image);

	//convert to B/W image
	unsigned char* data = image->getData();
	for (int i = 0; i < width*height; i++) if(data[i] != 0) data[i] = 255;
	
	//local blob counter to differentiate blobs
	unsigned char value = 254;

	// Find the Paper blob(s) in the copy (uses IntensityImage->integrate to get pos and size, overwrite values with a certain number) - see example in Bloblist.cc 138 and Blob.cc 24
	for (int i = 0; i < width*height; i++) if (data[i] == 255) try {

		tmpblobs->push_back( Blob( image, Point(i%width,i/width), value, id, minsize, maxsize ) );
/*			// integrate the spot and fill it with the current counter value
		size = (int)image->integrate( Point(i%width,i/width), pos, axis1, axis2, 255, value );
		
		// if the spot is too small or too big, wipe it out and abort
		if ( (size < minsize) || ((maxsize != 0) && (size > maxsize)) ){
			image->integrate( Point(i%width,i/width), pos, axis1, axis2, value, 0 );
			continue;
		}*/

		// adjust counter
		value--;
		id++;

		// did the frame-local blob counter overflow?
		if (value == 0) {
			value = 254;
			std::cerr << "Warning: too many blobs in AquaTopGenerator!" << std::endl;
		}
	} catch (...) { }

	// use the imagecopy to find all pixels outside of the blobs in the original image and overwrite them with 0
	for (int i = 0; i < width*height; i++)	if(data[i] == 0) tmpimagedata[i] = 0;

																							//copy shortimage to maskimage
																							//*maskimage = *shortimage;
																							//unsigned short* maskimagedata = maskimage->getSData();

	paperblobcounter = 0;
	hasnonpaperblob = false;
	bool ispaperblob = false;

	 //only keep the original values in the area of the blobs which are paper blobs			// (delete all paper blobs from mask image)
	for(int j = 254; j > value; j--)
	{
		ispaperblob = false;
		blobmin = 65535; blobmax = 0;
		for(int i = 0; i < width*height; i++) if(data[i] == j)
		{
			blobtmp = tmpimagedata[i];
			if(blobtmp < blobmin) blobmin = blobtmp;
			if(blobtmp > blobmax) blobmax = blobtmp;
		}
		//std::cout << "blobmin " << blobmin << " blobmax " << blobmax << " blobdiff " << blobmax - blobmin<< std::endl;
		//std::cout << "blobdiff " << blobmax - blobmin << " for value " << static_cast<unsigned>(value) << std::endl;
		if(blobmax - blobmin < paperdepthdiff) //--> potential paper blob
		{												
			// save the corresponding Blob in tmpblobs to paperblobs list
			for(std::vector<Blob>::iterator it = tmpblobs->begin(); it != tmpblobs->end(); it++)
			{
				if(it->value == j)
				{
					//std::cout << it->size << " =? " << it->axis1.length()*it->axis2.length()*4 << " diff: " << fabs((it->axis1.length() * it->axis2.length() * 4) - it->size) << std::endl;
					//<< " a1 length: " << it->axis1.length() << " a2 length: " << it->axis2.length() << std::endl;
					if(fabs((it->axis1.length() * it->axis2.length() * 4) - it->size) < 400) // --> paperblob
					{
						// use houghtransformation to check if the blob is a rectangle, cornerpoints are saved in papercorners if true is returned
						if(image->isRectangle(j, testimage1, testimage2, testimage3, testimage4, testimage5, houghTr, houghTtheta, houghTanglediff, it->axis2.length() / it->axis1.length(), papercorners))
						{

							//for (int i = 0; i < width*height; i++)	if(data[i] == j) maskimagedata[i] = 0;
							
							paperblobcounter++;
							ispaperblob = true;
							//std::cout << "paperblob??" << std::endl;
							
							paperblobs->push_back(Blob(*it));
							break;
						}
					}
				}
			}
			
		}
		if(!ispaperblob)//--> non-paper blob
		{
			hasnonpaperblob = true;
			for (int i = 0; i < width*height; i++)	if(data[i] == j) tmpimagedata[i] = 0;
		}
	}

	// Generate BG for the whole surface (Least-Squares Plane fitting) if at least one paper blob was found. Otherwise don't upgrade background image.
	if(paperblobcounter > 0)
	{
		GenerateBackground();
		*background = *tmpimage;
	}

	// BG subtraction
	*tmpimage = *(input->getShortImage());
	background->subtract( *tmpimage, *shortimage, invert );
	
	if(paperdepthdiff == 501) *shortimage = *testimage1;
	else if(paperdepthdiff == 502) *shortimage = *testimage2;
	/*else if(paperdepthdiff == 503) *shortimage = *testimage3;
	else if(paperdepthdiff == 504) *shortimage = *testimage4;
	else if(paperdepthdiff == 505) *shortimage = *testimage5;*/

	return 0;
}

// send blob list via OSC as TUIO 2.0
void AquaTopBGGenFilter::send( TUIOOutStream* oscOut ) {
	if(!hasnonpaperblob) oscOut->sendMessage("updatePapers");
	else oscOut->sendMessage("doNotUpdatePapers");

	std::vector<::Vector>::iterator cit = papercorners->begin();
	std::vector<::Vector>* corners = new std::vector<::Vector>;

	for (std::vector<Blob>::iterator it = paperblobs->begin(); it != paperblobs->end(); it++) {

		BasicBlob tmp = *it;
		tmp.type = 100; //paperblobs

	/*	std::cout << tmp.pos.x << " " << tmp.pos.y << " paperblob " << tmp.id << std::endl;
		if(cit != papercorners->end()) { std::cout << cit->x << " " << cit->y << " corner" << std::endl; cit++; }
		if(cit != papercorners->end()) { std::cout << cit->x << " " << cit->y << " corner" << std::endl; cit++; }
		if(cit != papercorners->end()) { std::cout << cit->x << " " << cit->y << " corner" << std::endl; cit++; }
		if(cit != papercorners->end()) { std::cout << cit->x << " " << cit->y << " corner" << std::endl; cit++; }*/

		tmp.pos.x  = tmp.pos.x  / (double)width; tmp.pos.y  = tmp.pos.y  / (double)height;
		tmp.peak.x = tmp.peak.x / (double)width; tmp.peak.y = tmp.peak.y / (double)height;

		/*if (hflip) {
			tmp.pos.x  = 1.0 - tmp.pos.x;
			tmp.peak.x = 1.0 - tmp.peak.x;
		}

		if (vflip) {*/
			tmp.pos.y  = 1.0 - tmp.pos.y;
			tmp.peak.y = 1.0 - tmp.peak.y;
		//}*/

		*oscOut << tmp;

		corners->clear();
		if(cit != papercorners->end()) 
		{ 
			cit->x /= (double) width; cit->y /= (double) height; cit->y = 1.0 - cit->y; corners->push_back(*cit); cit++;
			cit->x /= (double) width; cit->y /= (double) height; cit->y = 1.0 - cit->y; corners->push_back(*cit); cit++;
			cit->x /= (double) width; cit->y /= (double) height; cit->y = 1.0 - cit->y; corners->push_back(*cit); cit++;
			cit->x /= (double) width; cit->y /= (double) height; cit->y = 1.0 - cit->y; corners->push_back(*cit); cit++;
			oscOut->sendRectangle(tmp.id, *corners);
		}
	}
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
	case 5:
		OptionName = "Hough Tr";
		break;
	case 6:
		OptionName = "Hough Ttheta";
		break;
	case 7:
		OptionName = "Hough Tanglediff";
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
	case 5:
		OptionValue = houghTr;
		break;
	case 6:
		OptionValue = houghTtheta;
		break;
	case 7:
		OptionValue = houghTanglediff;
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
	case 5:
		if(overwrite) {
			houghTr = (delta < 0) ? 0 : (delta > 255) ? 255 : delta;
		} else {
			houghTr += delta;
			houghTr = (houghTr < 0) ? 0 : (houghTr > 255) ? 255 : houghTr;
		}
		break;
	case 6:
		if(overwrite) {
			houghTtheta = (delta < 0) ? 0 : (delta > width) ? width : delta;
		} else {
			houghTtheta += delta;
			houghTtheta = (houghTtheta < 0) ? 0 : (houghTtheta > width) ? width : houghTtheta;
		}
		break;
	case 7:
		if(overwrite) {
			houghTanglediff = (delta < 0) ? 0 : (delta > width) ? width : delta;
		} else {
			houghTanglediff += delta;
			houghTanglediff = (houghTanglediff < 0) ? 0 : (houghTanglediff > width) ? width : houghTanglediff;
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
	XMLNode->SetAttribute( "HoughTr",  houghTr);
	XMLNode->SetAttribute( "HoughTtheta",  houghTtheta);
	XMLNode->SetAttribute( "HoughTanglediff",  houghTanglediff);
	
	return XMLNode;
}

void AquaTopBGGenFilter::GenerateBackground()
{
	unsigned short* tmpimagedata = tmpimage->getSData();

	// calculate best fitting plane equation parameters: z = ax + by + c
	/* http://www.ahinson.com/algorithms_general/Sections/InterpolationRegression/RegressionPlanar.pdf
		(a,b,c) = m^-1 * v;
		m = 	sum(x^2)	sum(x*y)	sum(x)
				sum(x*y)	sum(y^2)	sum(y)
				sum(x)		sum(y)		n
		
		v = ( sum(x*z), sum(y*z), sum(z) )
	*/
	// 1) calculate m and v
	long long int sum_x = 0, sum_y = 0, sum_z = 0, sum_xx = 0, sum_yy = 0, sum_xy = 0, sum_xz = 0, sum_yz = 0;
	int z;
	int n = 0;
	for(int x = 0; x < width; x++)
	{
		for(int y = 0; y < height; y++)
		{
			if(tmpimagedata[y*width + x] != 0)
			{
				z = tmpimagedata[y*width + x];
				sum_x += x;	sum_y += y;	sum_z += z;	sum_xx += x*x; sum_yy += y*y; sum_xy += x*y; sum_xz += x*z;	sum_yz += y*z; n++;
			}
		}
	}
	Matrix m = { sum_xx, sum_xy, sum_x, sum_xy, sum_yy, sum_y, sum_x, sum_y, n };
	::Vector v;
	v.x = sum_xz; v.y = sum_yz; v.z = sum_z;

	// 2) calculate m_inv ( inverse of m )
	/*	http://ardoris.wordpress.com/2008/07/18/general-formula-for-the-inverse-of-a-3x3-matrix/ 
	formula:
		m^-1 =	factor *	m11m22-m12m12		m20m12-m10m22		m10m12-m20m11
							m12m02-m01m22		m00m22-m20m02		m20m01-m00m12
							m01m12-m11m02		m10m02-m00m12		m00m11-m10m01

		factor = 1 / det = 1 / ( m00(m11m22-m21m12) - m10(m01m22-m21m02) + m20(m01m12-m11m02) )
	*/
	double det = m.m00*(m.m11*m.m22-m.m21*m.m12) - m.m10*(m.m01*m.m22-m.m21*m.m02) + m.m20*(m.m01*m.m12-m.m11*m.m02);
	if(det == 0) return;
	double factor = 1 / det;
	Matrix m_inv = { factor, factor, factor, factor, factor, factor, factor, factor, factor };
	m_inv.m00 *= m.m11*m.m22-m.m12*m.m12;	
	m_inv.m10 *= (m.m20*m.m12-m.m10*m.m22);		
	m_inv.m20 *= m.m10*m.m12-m.m20*m.m11;
	m_inv.m01 *= m.m12*m.m02-m.m01*m.m22;		
	m_inv.m11 *= m.m00*m.m22-m.m20*m.m02;		
	m_inv.m21 *= m.m20*m.m01-m.m00*m.m12;
	m_inv.m02 *= m.m01*m.m12-m.m11*m.m02;		
	m_inv.m12 *= m.m10*m.m02-m.m00*m.m12;		
	m_inv.m22 *= m.m00*m.m11-m.m10*m.m01;

	// 3) calculate (a,b,c) = m^-1 * v;
	::Vector planeparams;
	planeparams.x = m_inv.m00 * v.x + m_inv.m10 * v.y + m_inv.m20 * v.z;
	planeparams.y = m_inv.m01 * v.x + m_inv.m11 * v.y + m_inv.m21 * v.z;
	planeparams.z = m_inv.m02 * v.x + m_inv.m12 * v.y + m_inv.m22 * v.z;

	// generate Background using the plane equation (z = ax + by + c) with calculated planeparams
	for(int x = 0; x < width; x++)	for(int y = 0; y < height; y++)
			if(tmpimagedata[y*width + x] == 0)	
			{tmpimagedata[y*width + x] = planeparams.x * x + planeparams.y * y + planeparams.z;
	}
}

	/*
	debugging...
	if(paperdepthdiff == 501) 
	{
		std::cout << "new frame -------------------------------------------------------------" << std::endl;
		std::cout << "Matrix m" << std::endl;
		std::cout << m.m00 << ", " << m.m10 << ", " << m.m20 << std::endl;
		std::cout << m.m01 << ", " << m.m11 << ", " << m.m21 << std::endl;
		std::cout << m.m02 << ", " << m.m12 << ", " << m.m22 << std::endl;
		std::cout << std::endl;
		std::cout << "Vector v: " << v.x << ", " << v.y << ", " << v.z << std::endl;
		std::cout << std::endl;
		std::cout << "det = " << det << " -- factor = " << factor << std::endl;
		std::cout << std::endl;
		std::cout << "Matrix m_inv" << std::endl;
		std::cout << m_inv.m00 << ", " << m_inv.m10 << ", " << m_inv.m20 << std::endl;
		std::cout << m_inv.m01 << ", " << m_inv.m11 << ", " << m_inv.m21 << std::endl;
		std::cout << m_inv.m02 << ", " << m_inv.m12 << ", " << m_inv.m22 << std::endl;
		std::cout << std::endl;
	}*/