/*************************************************************************\
*    Part of the TISCH framework - see http://tisch.sourceforge.net/      *
*  Copyright (c) 2006 - 2010 by Florian Echtler, TUM <echtler@in.tum.de>  *
*   Licensed under GNU Lesser General Public License (LGPL) 3 or later    *
\*************************************************************************/

#ifndef _TUIOSTREAM_H_
#define _TUIOSTREAM_H_

#include <string>
#include <vector>

#include <osc/OscOutboundPacketStream.h>
#include <ip/UdpSocket.h>
#include <BasicBlob.h>
#include <tisch.h>

#define TUIOSTREAM_BUFFER_SIZE 0x10000


class TUIOStream {

	template <class T> friend TUIOStream& operator<< ( TUIOStream& s, const T& t );

	public:

		TUIOStream( const char* target = "127.0.0.1", int port = TISCH_PORT_CALIB );

		void start();
		void setPrefix( std::string& _prefix );
		void send();

	protected:

		char buffer[TUIOSTREAM_BUFFER_SIZE];
		osc::OutboundPacketStream oscOut;
		UdpTransmitSocket transmitSocket;
		std::vector<int> alive;
		std::string prefix;
		int frame;
};


template <class T> TUIOStream& operator<< ( TUIOStream& s, const T& t ) {
	s.oscOut << t;
	return s;
}

template <> TUIOStream& operator<< <BasicBlob> ( TUIOStream& s, const BasicBlob& b );

#endif // _TUIOSTREAM_H_
