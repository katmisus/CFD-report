#ifndef _BOUND_COND_H_
#define _BORUND_COND_H_
#include "Types.h"
#include "Domain.h"

void BoundCond(Field& W, const Domain& dom);
void BoundCond(Field& W);
void ExchangeGhostCells(Field& W, const Domain& dom);

#endif
