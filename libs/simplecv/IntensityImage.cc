/*************************************************************************\
*    Part of the TISCH framework - see http://tisch.sourceforge.net/      *
*   Copyright (c) 2006 - 2011 by Florian Echtler <floe@butterbrot.org>    *
*   Licensed under GNU Lesser General Public License (LGPL) 3 or later    *
\*************************************************************************/

#include "IntensityImage.h"
#include "mmx.h"

#include <nanolibc.h>

#include <stdint.h>
#include <string.h>
#include <math.h>

#include <algorithm>
#include <fstream>
#include <string>

#include <Calibration.h>


/*inline int IntensityImage::pixelOffset(int x , int y, int channel) const { return (y*width)+x; }
inline unsigned char IntensityImage::getPixel(int x, int y, int channel) const { return data[(y*width)+x]; }
inline void IntensityImage::setPixel(int x, int y, unsigned char value, int channel) { data[(y*width)+x] = value; }*/


void IntensityImage::cross( int x, int y, int size, unsigned char color ) {
  for ( int i = x-(size/2); i <= x+(size/2); i++ ) setPixel(i,y,color);
  for ( int i = y-(size/2); i <= y+(size/2); i++ ) setPixel(x,i,color);
}

void IntensityImage::box( int x1, int y1, int x2, int y2, unsigned char color ) {
	for ( int i = x1; i <= x2; i++ ) { setPixel(i,y1,color); setPixel(i,y2,color); }
	for ( int i = y1; i <= y2; i++ ) { setPixel(x1,i,color); setPixel(x2,i,color); }
}


IntensityImage::IntensityImage( const char* path ) {

	int fwidth,fheight,fbpp;
	std::string magic,tmp;

	// open file with whitespace skipping
	std::ifstream myfile( path, std::ios::in );
	myfile >> std::skipws;

	// parse the header
	myfile >> magic;   myfile.ignore(1); if (myfile.peek() == '#') getline( myfile, tmp );
	myfile >> fwidth;  myfile.ignore(1); if (myfile.peek() == '#') getline( myfile, tmp );
	myfile >> fheight; myfile.ignore(1); if (myfile.peek() == '#') getline( myfile, tmp );
	myfile >> fbpp;

	if ((magic != "P5") || (fbpp > 255) || (fbpp < 1)) 
		throw std::runtime_error( std::string("IntensityImage: ") + std::string(path) + std::string(": no valid PGM file") );

	// init the base class
	init( fwidth, fheight, 1, 0, 0 );

	// skip one byte, read the rest
	myfile.ignore( 1 );
	myfile.read( (char*)data, size );
	myfile.close( );
}


int IntensityImage::histogram( int hg[] ) const {

	int max = 0;
	int pos = 0;

	for (int i = 0; i < 256; i++) hg[i] = 0;

	for (int i = 0; i < count; i++) hg[data[i]] += 1;

	for (int i = 0; i < 256; i++) 
		if (hg[i] > max) {
			max = hg[i];
			pos = i;
		}

	return pos;
}

int IntensityImage::intensity( ) const {
	#ifndef NOMMX
		return mmxintensity( data, count );
	#else
		long long int tmp = 0;
		for (int i = 0; i < count; i++) tmp += data[i];
		return (int)(tmp/(long long int)count);
	#endif
}


IntensityImage::IntensityImage( const IntensityImage& img, int downsample ):
	Image( (img.width)/2, (img.height)/2, 1 ) 
{
	int w1 = (img.width)/2;
	int h1 = (img.height)/2;
	for (int w = 0; w < w1; w++) for (int h = 0; h < h1; h++) {
		int tmp = 0;
		int offset = img.pixelOffset(w*2,h*2);
		tmp += img.data[offset  ];
		tmp += img.data[offset+1]; offset += img.width;
		tmp += img.data[offset  ];
		tmp += img.data[offset+1];
		setPixel(w,h,tmp/4,0);
	}
}


inline long long int   sum( int x ) { long long int t = x; return (t*(t+1))/2;         } // 1 + 2 + ... + x
inline long long int sqsum( int x ) { long long int t = x; return (t*(t+1)*(2*t+1))/6; } // 1^2 + 2^2 + ... + x^2

inline long long int   sum( int x1, int x2 ) { return   sum(x2)-  sum(x1-1); } // x1 + x1+1 + ... + x2
inline long long int sqsum( int x1, int x2 ) { return sqsum(x2)-sqsum(x1-1); } // x1^2 + (x1+1)^2 + ... + x2^2


void IntensityImage::scanSpan( int add, int dir, int x1, int x2, int y, Moments* m ) {

	// error condition:
	// scanspan: add -2, dir 1, x1 751, x2 752, y 436
	// scanspan: add 1, dir 1, x1 752, x2 751, y 437
	// scanspan: add 1, dir 1, x1 752, x2 751, y 437
	// ..... ad nauseam.
	// printf("scanspan: add %d, dir %d, x1 %d, x2 %d, y %d\n",add,dir,x1,x2,y);

	// add error checking for x1,x2 >= width
	if (x2 < x1) std::swap(x1,x2); // return
	if (x1 >= width) x1 = width-1;
	if (x2 >= width) x2 = width-1;

	if (add > 0) {

		long long int pcount = x2 - x1 + 1;
		long long int tmpsum = sum(x1,x2);
		int y0 = y-dir;

		m->m00 += pcount;
		m->m01 += pcount * y0;
		m->m02 += pcount * y0*y0;

		m->m10 += tmpsum;
		m->m20 += sqsum(x1,x2);
		m->m11 += tmpsum*y0;
	}

	if ((y < 0) || (y >= height)) return;

	register unsigned char oldcol = m->ov;
	register unsigned char newcol = m->nv;

	register int offs = width*y;
	register int x = x1;
	register int l = -1;

	// check for leaks to the left
	while ((x >= 0) && (data[offs+x] == oldcol)) { data[offs+x] = newcol; l = x--; }

	if (l != -1) {
		// leak found, see how far it extends to the right
		x = x1+1; while ((x < width) && (data[offs+x] == oldcol)) { data[offs+x] = newcol; x++; }
		scanSpan( -1, -dir, l, x-1, y-dir, m ); // FIXME: x-1 or x1-1?
	}

	// loop through normal spans
	while (x <= x2) {
		if (l != -1) scanSpan( 1, dir, l, x-1, y+dir, m );
		while ((x <= x2) && (data[offs+x] != oldcol)) x++;
		if (x <= x2) {
			l = x;
			while ((x < width) && (data[offs+x] == oldcol)) { data[offs+x] = newcol; x++; }
		} else {
			l = -1;
		}
	}

	// account for possible leak to the right
	if ((l != -1) && (x > x2)) { 
		scanSpan( -2, -dir, x2+1, x-1, y-dir, m );
		scanSpan(  2,  dir,    l, x-1, y+dir, m );
	}
}


long long int IntensityImage::integrate( Point start, Vector& centroid, Vector& axis1, Vector& axis2, unsigned char oldcol, unsigned char newcol, std::vector<Point>* border ) {

	Moments m = { 0, 0, 0, 0, 0, 0, oldcol, newcol };

	//m.m00 = integrate( start.x, start.y, &m );

	scanSpan( 0, -1, start.x, start.x, start.y, &m );
	scanSpan( 0,  1, start.x, start.x, start.y, &m );

	//std::cout << "moments: m00 = " << m00 << " m10 = " << m10 << " m01 = " << m01 << " m11 = " << m11 << " m20 = " << m20 << " m02 = " << m02 << std::endl;
	
	// calculate centroid: 1st order central moments
	double cx = (double)m.m10/(double)m.m00; centroid.x = cx;
	double cy = (double)m.m01/(double)m.m00; centroid.y = cy;

	// 2nd order central moments
	double mu20 = (double)m.m20/(double)m.m00 - cx*cx;
	double mu11 = (double)m.m11/(double)m.m00 - cx*cy;
	double mu02 = (double)m.m02/(double)m.m00 - cy*cy;

	// calculate eigenvalues lambda1,2 of image covariance matrix [ (mu20 mu11) (mu11 mu02) ]
	double tmp1 = mu20 + mu02;
	double tmp2 = mu20 - mu02;
	double tmp3 = sqrt( 4.0*mu11*mu11 + tmp2*tmp2 );

	double lambda1 = (tmp1 + tmp3)/2.0;
	double lambda2 = (tmp1 - tmp3)/2.0;

	// calculate eigenvectors using first matrix equation y = x * (lambda1,2 - mu20) / mu11
	axis1 = Vector( 1.0, ((lambda1-mu20)/mu11) );
	axis2 = Vector( 1.0, ((lambda2-mu20)/mu11) );

	// if result = NaN: repeat calculation using second equation y = x * mu11 / (lambda1,2 - mu02)
	if (isnan(axis1.y)) axis1.y = (mu11/(lambda1-mu02)); 
	if (isnan(axis2.y)) axis2.y = (mu11/(lambda2-mu02)); 

	// if result = Inf: normalize to base vector
	if (isinf(axis1.y)) axis1.set( 0.0, 1.0 );
	if (isinf(axis2.y)) axis2.set( 0.0, 1.0 );

	// normalize to axis lengths
	axis1.normalize(); axis1 = axis1 * sqrt( 3.0*fabs(lambda1) );
	axis2.normalize(); axis2 = axis2 * sqrt( 3.0*fabs(lambda2) );

	/* alternative angle calculation:
	double t = (mu02-mu20)/(2.0*mu11);
	double angle2 = atan(t+sqrt(t*t+1));*/
	centroid.z = 0.5 * atan( 2.0*mu11 / ( mu20 - mu02 ) );

	if (!border) return m.m00;

	// scan the border
	int neighbor[] = { -width-1, -width, -width+1, +1, width+1, width, width-1, -1 };

	int offs0 = start.y * width + start.x;
	int offset = offs0;
	int found = 0;

	while ((offset != offs0) || (!found)) {
		found = 0;
		for (int i = 0; i < 8; i++)
			if ((data[offset+neighbor[i]] == newcol) && (data[offset+neighbor[(i+1)%8]] != newcol)) {
				border->push_back( Point(offset%width,offset/width) );
				offset += neighbor[i];
				found = 1;
				break;
			}
		if (!found) { offset++; offs0 = offset; }
	}

	return m.m00;
}


void IntensityImage::despeckle( IntensityImage& target, unsigned char threshold ) const {
	#ifndef NOMMX
		memset( target.data, 0, width );
		memset( target.data+(width*(height-1)), 0, width );
		for (int j = 0; j < width-7; j+=6) mmxdespeckle( data+j, target.data+j, height, width, threshold );
		mmxdespeckle( data+width-7, target.data+width-7, height, width, threshold );
	#else
		int neighbor[] = { -width-1, -width, -width+1, -1, 0, +1, width-1, width, width+1 };
		int sum;
		target.clear();
		for (int i = width+1; i < count-width-1; i++) {
			sum = 0;
			for (int j = 0; j < 9; j++) if (data[i+neighbor[j]]) sum++;
			target.data[i] = (sum >= threshold) ? 255 : 0;
		}
	#endif
}

void IntensityImage::lowpass( IntensityImage& target, unsigned char range, unsigned char mode ) const {
	int sum;
	target.clear();
	if(mode == 0) //horizontal
	{
		for (int i = range; i < count-range; i++) 
		{
			sum = 0;
			for (int j = -range; j <= range; j++) if (data[i+j]) sum++;
			target.data[i] = (sum == 2*range + 1) ? 255 : 0;
		}
	}
	else if (mode == 1) //vertical
	{
		for (int i = range*width; i < count-range*width; i++) 
		{
			sum = 0;
			for (int j = -range; j <= range; j++) if (data[i+j*width]) sum++;
			target.data[i] = (sum == 2*range + 1) ? 255 : 0;
		}
	}
	else if (mode == 2) //horizontal + vertical
	{
		for (int i = range*width + range; i < count-range*width - range; i++) 
		{
			sum = 0;
			for (int j = -range; j <= range; j++) {
				if (data[i+j]) sum++; 
				if (data[i+j*width]) sum++;
			}
			target.data[i] = (sum == 4*range + 2) ? 255 : 0;
		}
	}
}


void IntensityImage::houghLine( IntensityImage& target ) const {

	float theta,rho;
	unsigned char tmp;
	int x, y, j, offset, res, max = 0;

	int  size = target.width * target.height;
	int* accu = new int[ size ];
	for (int i = 0; i < size; i++) accu[i] = 0;

	float dt = M_PI/(float)target.width;

	float rf = 2.0*sqrtf(width*width+height*height);

	//float cos_t;
	//float sin_t;

	for (x = 0; x < width; x++) for (y = 0; y < height; y++) {
	
		tmp = data[ x+y*width ];
		if (tmp == 0) continue;

		for (j = 0, theta = 0; j < target.width; j++,theta+=dt) {
			//cos_t = cosf(theta);
			//sin_t = sinf(theta);
			//if(abs(cos_t) < 0.1) {cos_t = 0; sin_t = 1;}
			//if(abs(sin_t) < 0.1) {cos_t = 1; sin_t = 0;}
			rho = x * cosf(theta) + y * sinf(theta);//x * cos_t + y * sin_t;//
			rho = (rho/rf+0.5)*target.height;

			offset = j+target.width*(int)roundf(rho);

			res = accu [ offset ] + 1;
			if (res > max) max = res;
			accu [ offset ] = res;
		}
	}

	for (int i = 0; i < size; i++) target.data[i] = (unsigned char) roundf( ((float)accu[i] * 255.0) / (float)max );

	delete[] accu;
}

void IntensityImage::areamask( IntensityImage& target, std::vector<int> edgepoints) const
{
	if( edgepoints.empty() )
		memcpy(target.data,data,size);
	else
		memset(target.data, 0, size);
		for(std::vector<int>::iterator it = edgepoints.begin(); it != edgepoints.end(); it = it + 2)
			memcpy((target.data +  (*it)), (data + (*it)), (*(it+1) - *it)); 
}

void IntensityImage::gradient( char* xgrad, char* ygrad ) {

	int xval,yval,offs;

	for (int x = 1; x < width-1; x++) for (int y = 1; y < height-1; y++) {
		
		offs = pixelOffset(x,y);
		
		xval = getPixel(x+1,y) - getPixel(x-1,y);
		yval = getPixel(x,y-1) - getPixel(x,y+1);
		
		xgrad[offs] = xval/2;
		ygrad[offs] = yval/2;
	}
}

void IntensityImage::mask( IntensityImage& target ) {
	for (int i = 0; i < size; i++) if (data[i] == 0) target.data[i] = 0;
}
	
void IntensityImage::sobel( unsigned char* target ) {

	int xval,yval,value;

	for (int x = 1; x < width-1; x++) for (int y = 1; y < height-1; y++) {

		xval =   getPixel(x+1,y-1) -   getPixel(x-1,y-1)
				 + 2*getPixel(x+1,y  ) - 2*getPixel(x-1,y  )
				 +   getPixel(x+1,y+1) -   getPixel(x-1,y+1);

		yval =   getPixel(x-1,y-1) -   getPixel(x-1,y+1)
				 + 2*getPixel(x  ,y-1) - 2*getPixel(x  ,y+1)
				 +   getPixel(x+1,y-1) -   getPixel(x+1,y+1);

		value = (int)sqrt((double)(xval*xval+yval*yval));
		value /= 4;
		value = (value > 255) ? 255 : value;
		value = (value < 0) ? 255 : value;
		target[pixelOffset(x,y)] = value;
	}
}

void IntensityImage::sobel( IntensityImage& target ) {
	sobel( target.data );
}

void IntensityImage::sobel() {
	// evil pointer shuffling
	unsigned char* tmp = new unsigned char[size];
	sobel( tmp );
	delete[] data;
	data = tmp;
}


void IntensityImage::adaptive_threshold( int radius, int bias, IntensityImage& target ) const {

	int* colsum = new int[width];
	int colcnt = radius+1;

	// initialize column sums
	for (int x = 0; x < width; x++) {
		colsum[x] = 0;
		for (int y = 0; y <= radius; y++)
			colsum[x] += getPixel(x,y);
	}

	for (int y = 0; y < height; y++) {

		for (int x = 0; x < width; x++) {

			// add the column sums from -radius to +radius
			int thr = 0, cnt = 0;
			for (int i = x-radius; i <= x+radius; i++) {
				if ((i < 0) || (i >= width)) continue;
				thr += colsum[i];
				cnt += colcnt;
			}

			// calculate threshold from mean and bias
			thr /= cnt;
			thr -= bias;

			// apply to pixel
			target.setPixel( x, y, ( (getPixel(x,y) > thr) ? 255 : 0 ) );
		}

		// update column sums
		int lower = y-radius-1;
		int upper = y+radius+1;

		if (lower >=    0) { colcnt--; for (int x = 0; x < width; x++) colsum[x] -= getPixel( x, lower ); }
		if (upper < width) { colcnt++; for (int x = 0; x < width; x++) colsum[x] += getPixel( x, upper ); }
	}
	delete colsum;
}

int IntensityImage::threshold( unsigned char value, unsigned char minvalue ) {
	return threshold( value, *this, minvalue );
}

int IntensityImage::threshold( unsigned char value, IntensityImage& target , unsigned char minvalue ) const {
	#ifndef NOMMX
		mmxthreshold( data, target.data, size, value, minvalue );
		return 0;
	#else
		unsigned char tmp;
		int hits = 0;
		for (int i = 0; i < count; i++) {
			tmp = ((data[i] > value) && (data[i] <= minvalue)) ? 255 : 0;
			if (tmp) hits++;
			target.data[i] = tmp;
		}
		return hits;
	#endif
}


void IntensityImage::undistort( Vector scale, Vector delta, double coeff[5], IntensityImage& target ) const {

	Vector temp;
	target.clear();

	for (int u = 0; u < width; u++) for (int v = 0; v < height; v++) {

		temp = Vector( u, v, 0 );
		::undistort( temp, scale, delta, coeff );

		if ((temp.x < 0) || (temp.x >=  width)) continue;
		if ((temp.y < 0) || (temp.y >= height)) continue;
		target.setPixel( u, v, getPixel( (int)temp.x, (int)temp.y ) );
	}
}


void IntensityImage::invert( ) {
	invert( *this );
}

void IntensityImage::invert( IntensityImage& target ) const {
	unsigned char tmp;
	for (int i = 0; i < count; i++) {
		tmp = 255-data[i];
		target.data[i] = tmp;
  }
}


IntensityImage& IntensityImage::operator-=( const IntensityImage& arg ) {
	int tmp;
	for (int i = 0; i < count; i++) {
		tmp = (int)data[i] - (int)arg.data[i];
		data[i] = (tmp > 0) ? ((tmp < 255) ? tmp : 255 ) : 0;
	}
	return *this;
}


void IntensityImage::subtract( const IntensityImage& i1, const IntensityImage& i2 ) {
	int tmp;
	for (int i = 0; i < count; i++) {
		tmp = (int)i1.data[i] - (int)i2.data[i];
		data[i] = (tmp > 0) ? ((tmp < 255) ? tmp : 255 ) : 0;
	}
}


std::ostream& operator<<( std::ostream& s, const IntensityImage& i ) {
	s << "P5 " << i.width << " " << i.height << " 255 ";
	s.write( (char*)i.data, i.size );
	return s;	
}


// adapted from EquisFtir/touch/filter_band_mid.cpp

#define clamp(v,l,u) (v<l?l:(v>u?u:v))

struct BandMid {
	int outer;
	int inner;

	int rows;
	int cols;

	int shift;
	int im;

	void apply_row(const uint8_t* src, uint8_t* dst);
	void apply_col(const uint8_t* src, uint8_t* dst);
};

void BandMid::apply_row(const uint8_t* src, uint8_t* dst) {
	
	long a0 = 0;
	long a1 = 0;
	long a2 = 0;

	for (int c = -outer; c < 0; ++c) {
		a0 += a1; a1 += a2;

		if (c >= -outer) a2 -= 1L * src[c + outer];
		if (c >= -outer/2) a2 += 2L * src[c + outer/2];

		if (c >= -inner) a2 += 1L * im * src[c + inner];
		if (c >= -inner/2) a2 -= 2L * im * src[c + inner/2];
	}

	for (int c = 0; c < outer; ++c) {
		dst[c] = static_cast<uint8_t>(clamp(a0 >> shift, 0L, 255L));

		a0 += a1; a1 += a2;

		if (c < cols - outer) a2 -= 1L * src[c + outer];
		if (c < cols - outer/2) a2 += 2L * src[c + outer/2];
		if (c >= outer/2) a2 -= 2L * src[c - outer/2];
		if (c >= outer) a2 += 1L * src[c - outer];

		if (c < cols - inner) a2 += 1L * im * src[c + inner];
		if (c < cols - inner/2) a2 -= 2L * im * src[c + inner/2];
		if (c >= inner/2) a2 += 2L * im * src[c - inner/2];
		if (c >= inner) a2 -= 1L * im * src[c - inner];
	}

	for (int c = outer; c < cols - outer; ++c) {
		dst[c] = static_cast<uint8_t>(clamp(a0 >> shift, 0L, 255L));

		a0 += a1; a1 += a2;

		a2 -= 1L * src[c + outer];
		a2 += 2L * src[c + outer/2];
		a2 -= 2L * src[c - outer/2];
		a2 += 1L * src[c - outer];

		a2 += 1L * im * src[c + inner];
		a2 -= 2L * im * src[c + inner/2];
		a2 += 2L * im * src[c - inner/2];
		a2 -= 1L * im * src[c - inner];
	}

	for (int c = cols - outer; c < cols; ++c) {
		dst[c] = static_cast<uint8_t>(clamp(a0 >> shift, 0L, 255L));

		a0 += a1; a1 += a2;

		if (c < cols - outer) a2 -= 1L * src[c + outer];
		if (c < cols - outer/2) a2 += 2L * src[c + outer/2];
		if (c >= outer/2) a2 -= 2L * src[c - outer/2];
		if (c >= outer) a2 += 1L * src[c - outer];

		if (c < cols - inner) a2 += 1L * im * src[c + inner];
		if (c < cols - inner/2) a2 -= 2L * im * src[c + inner/2];
		if (c >= inner/2) a2 += 2L * im * src[c - inner/2];
		if (c >= inner) a2 -= 1L * im * src[c - inner];
	}
}

void BandMid::apply_col(const uint8_t* src, uint8_t* dst) {

	long a0 = 0;
	long a1 = 0;
	long a2 = 0;

	for (int r = -outer; r < 0; ++r) {
		a0 += a1; a1 += a2;

		if (r >= -outer) a2 -= 1L * src[(r + outer) * cols];
		if (r >= -outer/2) a2 += 2L * src[(r + outer/2) * cols];

		if (r >= -inner) a2 += 1L*im * src[(r + inner) * cols];
		if (r >= -inner/2) a2 -= 2L*im * src[(r + inner/2) * cols];
	}

	for (int r = 0; r < outer; ++r) {
		dst[r*cols] = static_cast<uint8_t>(clamp(a0 >> shift, 0L, 255L));

		a0 += a1; a1 += a2;

		if (r < rows - outer) a2 -= 1L * src[(r + outer) * cols];
		if (r < rows - outer/2) a2 += 2L * src[(r + outer/2) * cols];
		if (r >= outer/2) a2 -= 2L * src[(r - outer/2) * cols];
		if (r >= outer) a2 += 1L * src[(r - outer) * cols];

		if (r < rows - inner) a2 += 1L*im * src[(r + inner) * cols];
		if (r < rows - inner/2) a2 -= 2L*im * src[(r + inner/2) * cols];
		if (r >= inner/2) a2 += 2L*im * src[(r - inner/2) * cols];
		if (r >= inner) a2 -= 1L*im * src[(r - inner) * cols];
	}

	for (int r = outer; r < rows - outer; ++r) {
		dst[r*cols] = static_cast<uint8_t>(clamp(a0 >> shift, 0L, 255L));

		a0 += a1; a1 += a2;

		a2 -= 1L * src[(r + outer) * cols];
		a2 += 2L * src[(r + outer/2) * cols];
		a2 -= 2L * src[(r - outer/2) * cols];
		a2 += 1L * src[(r - outer) * cols];

		a2 += 1L * im * src[(r + inner) * cols];
		a2 -= 2L * im * src[(r + inner/2) * cols];
		a2 += 2L * im * src[(r - inner/2) * cols];
		a2 -= 1L * im * src[(r - inner) * cols];
	}

	for (int r = rows - outer; r < rows; ++r) {
		dst[r*cols] = static_cast<uint8_t>(clamp(a0 >> shift, 0L, 255L));

		a0 += a1; a1 += a2;

		if (r < rows - outer) a2 -= 1L * src[(r + outer) * cols];
		if (r < rows - outer/2) a2 += 2L * src[(r + outer/2) * cols];
		if (r >= outer/2) a2 -= 2L * src[(r - outer/2) * cols];
		if (r >= outer) a2 += 1L * src[(r - outer) * cols];

		if (r < rows - inner) a2 += 1L*im * src[(r + inner) * cols];
		if (r < rows - inner/2) a2 -= 2L*im * src[(r + inner/2) * cols];
		if (r >= inner/2) a2 += 2L*im * src[(r - inner/2) * cols];
		if (r >= inner) a2 -= 1L*im * src[(r - inner) * cols];
	}
}


void IntensityImage::bandpass( IntensityImage& target, int outer, int inner ) const {

	BandMid s;
	s.outer = outer;
	s.inner = inner;

	s.rows = height;
	s.cols = width;

	// outer = 16, inner = 8, shift = 3, add 1
	// outer = 16, inner = 4, shift = 6, add 4
	s.im = 8;
	s.shift = static_cast<int>(log(0.25 * s.inner * s.inner * s.inner) / log(2.0)) + 1;

	uint8_t* tmp = new uint8_t[size];

	for (int r = 0; r < s.rows; r++) 
		s.apply_row( data+(r*width), tmp+(r*width) );

	for (int c = 0; c < s.cols; c++)
		s.apply_col( tmp+c, target.data+c );

	for (int r = 0; r < s.outer; r++) {
		for (int c = 0; c < s.outer; c++) {
			target.setPixel(c, r, 0);
			target.setPixel(s.cols - 1 - c, r, 0);
			target.setPixel(c, s.rows - 1 - r, 0);
			target.setPixel(s.cols - 1 - c, s.rows - 1 - r, 0);
		}
	}

	delete[] tmp;
}


bool IntensityImage::isRectangle( unsigned char value, IntensityImage* t1, IntensityImage* t2, IntensityImage* t3, IntensityImage* t4, IntensityImage* t5, int Tr , int Ttheta, int Tanglediff, double axisratio, std::vector<::Vector>* papercorners) const
{
	bool debugcout = false;
	// copy image
	//IntensityImage* tmp = new IntensityImage(width, height);
	*t1 = *this;//*tmp = *this;
	unsigned char* copieddata = t1->getData();//tmp->getData();

	// create B/W image where value->W (255), !value->B (0)
	for(int i = 0; i < width * height; i++)
	{
		copieddata[i] = (copieddata[i] == value) ? 255 : 0;
	}

	// sobel
	//IntensityImage* sobelimage = new IntensityImage(width, height);
	t1->sobel (*t2);//tmp->sobel( *sobelimage );

	// houghtransformation (use flo's houghline?)
	t2->houghLine(*t3);//sobelimage->houghLine( *sobelimage );

	t3->butterflyevaluator(*t4, 2, 2);

	// extract local maxima
	t4->localextrema(*t5,9,5);
	
	// search for rectangle in HT image (using algorithm of paper) (accumulated values in HT > Tc, Tc should be axisratio * 255)
	// find peaks
	unsigned char* enhancedhoughdata = t5->getData();
	int houghwidth = t3->getWidth();

	std::vector<int> tmppeaks;
	for(int i = 0; i < t5->getCount(); i++)	if(enhancedhoughdata[i] != 0) tmppeaks.push_back(i); 
	if(debugcout) std::cout << tmppeaks.size() << " tmppeaks in HT" << std::endl;
	
	// store all peaks were houghdata < Tc
	unsigned char* houghdata = t3->getData();
	int Tc = 100;//(int)(0.5*axisratio * 255);
	std::vector<int> peaks;
	for(int i = 0; i < tmppeaks.size(); i++) if(houghdata[tmppeaks[i]] > Tc) {peaks.push_back(tmppeaks[i]); t2->cross(tmppeaks[i]%houghwidth, (int)(tmppeaks[i]/houghwidth),10,255); if(debugcout) std::cout << tmppeaks[i]%houghwidth << " " << tmppeaks[i]/houghwidth << std::endl; }
	if(debugcout) std::cout << peaks.size() << " peaks in HT" << std::endl;
	if(peaks.size() > 6) return false;

	// find pairs with approximately same angle: delta theta = |theta1 - theta2| < Ttheta
	struct Pair { int p1; int p2; };
	std::vector<Pair> pairs;
	for(int i = 0; i < peaks.size() - 1; i++) for(int j = i+1; j < peaks.size(); j++)
	{
		int diff = abs(peaks[i]%houghwidth - peaks[j]%houghwidth);
		//diff = (diff > 0.5*houghwidth) ? houghwidth - diff : diff;
		int ydiff = abs(peaks[i]/houghwidth - peaks[j]/houghwidth);
		//if(ydiff > 15) ydiff = abs(peaks[i]/houghwidth - (height - peaks[j]/houghwidth));
		if( diff < Ttheta &&  ydiff > Tr)
		{
			Pair tmppair; tmppair.p1 = peaks[i]; tmppair.p2 = peaks[j];
			pairs.push_back(tmppair);
		}
	}

	if(debugcout) std::cout << pairs.size() << " pairs in HT" << std::endl;
	if(pairs.size() == 0) return false;

	int counter = 0;
	// find rectangles: alpha = 0.5 * (theta1 + theta2); delta alpha = ||alpha1 - alpha2| - 90°| < Talpha
	for(int i = 0; i < pairs.size() - 1; i++) for(int j = i+1; j < pairs.size(); j++)
	{
		double pairdist = pairs[i].p1%houghwidth - pairs[j].p1%houghwidth;//abs(0.5*(pairs[i].p1%width + pairs[i].p2%width) - 0.5*(pairs[j].p1%width + pairs[j].p2%width));
		pairdist = (pairdist > 0.5*houghwidth) ? houghwidth - pairdist : pairdist;
		double zeroanglediff = abs(pairdist - 0.5*houghwidth);
		if(zeroanglediff < Tanglediff) //tolerance for 90° difference between lines
		{
			counter++;
			//just use the first found rectangle... (this could be improved by considering all rectangles and calculate the estimated error of each and take the lowest)
			//calculate cornerpoints of the rectangle and save them (to send them to client)
			// pair[i].p1+2 = Indices of 2 houghdata points
			// pair[j].p1+2 = Indices of the other 2 houghdata points
			// calculate r and theta for these 4 points (theta = %houghwidth * angledistribution, r = /houghwidth)
			// these are the parameters r and theta in r = x * cos(theta) + y * sin(theta) for the 4 lines
			// use pairs of equations (one from pair[i], one from pair[j]) in a LGS to calculate each corner point
			int r1, r2, r3, r4;
			double theta1, theta2, theta3, theta4;
			r1 = pairs[i].p1 / houghwidth;
			r2 = pairs[i].p2 / houghwidth;
			r3 = pairs[j].p1 / houghwidth;
			r4 = pairs[j].p2 / houghwidth;
			//std::cout << r1 << " " << r2 << " " << r3 << " " << r4 << std::endl;
			theta1 = (pairs[i].p1 % houghwidth) * M_PI / houghwidth;
			theta2 = (pairs[i].p2 % houghwidth) * M_PI / houghwidth;
			theta3 = (pairs[j].p1 % houghwidth) * M_PI / houghwidth;
			theta4 = (pairs[j].p2 % houghwidth) * M_PI / houghwidth;
			//std::cout<< std::endl;
			/*::Vector cornerrr = calculateCornerPoint(r1,r3,theta1,theta3);
			papercorners->push_back(cornerrr);
			papercorners->push_back(cornerrr);
			papercorners->push_back(cornerrr);
			papercorners->push_back(cornerrr);*/
			papercorners->push_back(calculateCornerPoint(r1,r3,theta1,theta3));
			papercorners->push_back(calculateCornerPoint(r1,r4,theta1,theta4));
			papercorners->push_back(calculateCornerPoint(r2,r3,theta2,theta3));
			papercorners->push_back(calculateCornerPoint(r2,r4,theta2,theta4));
			for(std::vector<::Vector>::iterator it = papercorners->begin(); it != papercorners->end(); it++) t2->cross(it->x, it->y, 20, 255);
			return true;
		}
	}
	if(debugcout) std::cout << counter << " rectangles" << std::endl;
	//std::cout << "false" << std::endl;
	if(counter == 0) return false;
	return true;
}


void IntensityImage::butterflyevaluator( unsigned char* target, int w, int h, int targetsize) {

	int value, divisor;
	int* values= new int[ targetsize ];
	int max = 0;

	int tmpx, tmpy;

	//for (int x = w; x < width - w; x++) for (int y = h; y < height - h; y++) {
	for (int x = 0; x < width; x++) for (int y = h; y < height - h; y++) {
		value = 0;
		divisor = 0;
		for (int i = -w; i < w; i++) for (int j = -h; j < h; j++)
		{
			tmpx = x+i;
			tmpy = y+i;
			if(tmpx < 0) { tmpx = width + tmpx; tmpy = height - tmpy; }
			else if (tmpx >= width) {tmpx = tmpx - width; tmpy = height - tmpy; }
			divisor += getPixel(tmpx,tmpy);
		}
		if(divisor == 0) continue;
		value = (4 * h * w + 2) * getPixel(x,y) * getPixel(x,y);
		value /= divisor;
		if (value > max) max = value;
		values[pixelOffset(x,y)] = value;
	}

	for (int i = 0; i < size; i++) target[i] = (unsigned char) roundf( ((float)values[i] * 255.0) / (float)max );

	delete[] values;
}

void IntensityImage::butterflyevaluator( IntensityImage& target, int w, int h ) {
	butterflyevaluator( target.data , w, h, target.width * target.height);
}

void IntensityImage::butterflyevaluator( int w, int h ) {
	// evil pointer shuffling
	unsigned char* tmp = new unsigned char[size];
	butterflyevaluator( tmp, w, h, width * height );
	delete[] data;
	data = tmp;
}

void IntensityImage::localextrema( IntensityImage& target, int w, int h) const {

	for(int x = 0; x < target.width; x++) for(int y = 0; y < target.height; y++) target.setPixel(x,y,1);

	unsigned char tmp;
	bool ismax;
	int xpos, ypos;

	for (int x = 0; x < width; x++) for (int y = h; y < height - h; y++) 
	{
		if(target.getPixel(x,y) == 0) continue;
		tmp = data[ x+y*width ];
		if (tmp == 0) continue;
		ismax = true;
		for (int i = -w; i < w; i++) for (int j = -h; j < h; j++)
		{
			xpos = x+i;
			//(xpos < 0) ? ypos = height - (y+j) : ypos = (y+j+height)%height;
			//(xpos > width) ? ypos = height - (y+j) : ypos = (y+j+height)%height;
			if(xpos < 0 || xpos > width) ypos = height - (y+j);
			else ypos = y+j;//ypos = (y+j+height)%height;

			if(xpos < 0) xpos += width;
			else if(xpos >= width) xpos -= width;

			if(tmp < getPixel(xpos,ypos))
			{
				target.setPixel(x,y,0);
				//there is a higher value at xpos,ypos --> start recursive call
				if(target.getPixel(xpos,ypos) == 1) localextremaRecursive(target,w,h,xpos,ypos);				
				ismax = false;
				break;
			}
			else
			{
				target.setPixel(xpos,ypos,0);
			}
		}
		if(ismax) target.setPixel(x,y,tmp);
	}
}

void IntensityImage::localextremaRecursive( IntensityImage& target, int w, int h, int x, int y) const {

	unsigned char tmp;
	bool ismax;
	int xpos, ypos;

	tmp = data[ x+y*width ];
	if (tmp == 0 || target.getPixel(x,y) == 0) return;
	ismax = true;
	for (int i = 0; i < w; i++) for (int j = 0; j < h; j++)
	{
		xpos = x+i;
		if(xpos > width) ypos = height - (y+j);
		else ypos = y+j;

		if(xpos >= width) xpos -= width;

		if(tmp < getPixel(xpos,ypos))
		{
			target.setPixel(x,y,0);
			//there is a higher value at xpos,ypos --> start recursive call
			if(target.getPixel(xpos,ypos) == 1) localextremaRecursive(target,w,h,xpos,ypos);				
			ismax = false;
			break;
		}
		else
		{
			target.setPixel(xpos,ypos,0);
		}
	}
	if(ismax) target.setPixel(x,y,tmp);
}

// r = x * cos(theta) + y * sin(theta)
::Vector IntensityImage::calculateCornerPoint( int r1, int r2, double theta1, double theta2 ) const
{
	//std::cout << r1 << " " << r2 << " " << height << " " << theta1 << " " << theta2 << std::endl;
	float rf = 2.0*sqrtf(width*width+height*height); //this width and height is from the original image

	r1 = ( (double)r1 /  (0.5*height) - 0.5) * rf; //only works if height of hough-image == height of original image (when using r1 / height - 0.5). the height in this calculation should be the height of the hough image
	r2 = ( (double)r2 /  (0.5*height) - 0.5) * rf; //only works if height of hough-image == height of original image (when using r2 / height - 0.5). the height in this calculation should be the height of the hough image
	//std::cout << r1 << " " << r2 << " " << theta1 << " " << theta2 << std::endl;

	double cos_t1 = cos(theta1);
	double cos_t2 = cos(theta2);

	::Vector point;
	int x,y;
	x=1;y=1;
	if(abs(cos_t1) < 0.00001/*< 0.1*/) //theta = 90°, theta2 = 0° oder 180°
	{
		y = r1; //r1 / sin(theta1);
		x = abs(r2); //(r2 - y * sin(theta2)) / cos(theta2);
		//std::cout << "t1  " << x << " " << y << std::endl;
	}
	else if (abs(cos_t2) < 0.00001/*< 0.1*/)//theta2 = 90°, theta1 = 0° oder 180°
	{
		y = r2; //r2 / sin(theta2);
		x = abs(r1); //(r1 - y * sin(theta1)) / cos(theta1);
		//std::cout << "t2  " << x << " " << y << std::endl;
	}
	else
	{
		double par1 = r1 / cos(theta1);
		double par2 = r2 / cos(theta2);
		y = (par1 - par2) / (tan(theta1) - tan(theta2));
		x = par1 - y * tan(theta1);
	}
	//std::cout << "Corner: (" << x << ", " << y << ")" << std::endl;
	point.x = x;
	point.y = y;
	return point;
}
