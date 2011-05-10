/*************************************************************************\
*    Part of the TISCH framework - see http://tisch.sourceforge.net/      *
*  Copyright (c) 2006 - 2009 by Florian Echtler, TUM <echtler@in.tum.de>  *
*   Licensed under GNU Lesser General Public License (LGPL) 3 or later    *
\*************************************************************************/

#ifndef _MULTIBLOBROTATION_H_
#define _MULTIBLOBROTATION_H_

#include "Rotation.h"

class TISCH_SHARED MultiBlobRotation: public Rotation {

	public:

		 MultiBlobRotation( int tf = (1<<INPUT_TYPE_COUNT)-1 );
		~MultiBlobRotation();

		void load( InputState& state );

		const char* name() const { return "MultiBlobRotation"; }

};

#endif // _MULTIBLOBROTATION_H_
