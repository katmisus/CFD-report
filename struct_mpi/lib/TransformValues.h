#ifndef _TRANSFORM_VALUES_H_
#define _TRANSFORM_VALUES_H_
#include "Types.h"

void ConvertWtoU(Field W, Field& U, int dir);
void ConvertWtoU(State W, State& U);
void ConvertUtoW(Field& W, Field U, int dir);

#endif