/**
 *  HEVC Encoder
 *  - Marko Viitanen ( fador at iki.fi ), Tampere University of Technology, Department of Pervasive Computing.
 */

/*! \file search.c
    \brief searching
    \author Marko Viitanen
    \date 2013-04
    
    Search related functions
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "global.h"
#include "config.h"
#include "bitstream.h"
#include "picture.h"
#include "cabac.h"
#include "encoder.h"
#include "filter.h"
