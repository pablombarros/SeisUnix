/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

#ifndef SU_XDR
#define SU_XDR

#include <stdio.h>

#ifdef SUXDR
#include <rpc/types.h>
#include <rpc/xdr.h>
#endif


#include "su.h"
#include "segy.h"

#ifdef SUXDR
int xdrhdrsub(XDR *segyxdr, segy *tp);
int xdrbhdrsub(XDR *segyxdr, bhed *bhp);
#else
void xdrhdrsub();
void xdrbhdrsub();
#endif

#endif /* SU_XDR */ 

