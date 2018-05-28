/*
 *  Copyright (C) 2004-2017 Edward F. Valeev
 *
 *  This file is part of Libint.
 *
 *  Libint is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Libint.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <libint2.h>
#include <_elecpot_P_D_prereq.h>

#ifdef __cplusplus
extern "C" {
#endif
void _elecpot_P_D(const Libint_t* inteval) {

_libint2_static_api_bzero_short_(inteval->stack+0*1,19)
const int contrdepth = inteval->contrdepth;
#ifdef __INTEL_COMPILER
#pragma ivdep
#endif
for(int c = 0; c<contrdepth; c++) {
_elecpot_P_D_prereq(inteval+c, inteval->stack);
}

LIBINT2_REALTYPE*const  stack = inteval->stack;
{
const int hsi = 0;
{
const int lsi = 0;
{
const int vi = 0;
LIBINT2_REALTYPE fp33;
fp33 = inteval->AB_x[vi] * inteval->stack[((hsi*1+4)*1+lsi)*1];
LIBINT2_REALTYPE fp32;
fp32 = inteval->stack[((hsi*1+3)*1+lsi)*1] + fp33;
LIBINT2_REALTYPE fp35;
fp35 = inteval->AB_x[vi] * inteval->stack[((hsi*1+2)*1+lsi)*1];
LIBINT2_REALTYPE fp34;
fp34 = inteval->stack[((hsi*1+1)*1+lsi)*1] + fp35;
LIBINT2_REALTYPE fp81;
fp81 = inteval->AB_y[vi] * fp34;
LIBINT2_REALTYPE fp80;
fp80 = fp32 + fp81;
inteval->stack[((hsi*18+20)*1+lsi)*1] = fp80;
LIBINT2_REALTYPE fp29;
fp29 = inteval->AB_x[vi] * inteval->stack[((hsi*1+6)*1+lsi)*1];
LIBINT2_REALTYPE fp28;
fp28 = inteval->stack[((hsi*1+5)*1+lsi)*1] + fp29;
LIBINT2_REALTYPE fp79;
fp79 = inteval->AB_z[vi] * fp34;
LIBINT2_REALTYPE fp78;
fp78 = fp28 + fp79;
inteval->stack[((hsi*18+21)*1+lsi)*1] = fp78;
LIBINT2_REALTYPE fp31;
fp31 = inteval->AB_y[vi] * inteval->stack[((hsi*1+2)*1+lsi)*1];
LIBINT2_REALTYPE fp30;
fp30 = inteval->stack[((hsi*1+4)*1+lsi)*1] + fp31;
LIBINT2_REALTYPE fp77;
fp77 = inteval->AB_y[vi] * fp30;
LIBINT2_REALTYPE fp75;
fp75 = inteval->AB_y[vi] * inteval->stack[((hsi*1+4)*1+lsi)*1];
LIBINT2_REALTYPE fp74;
fp74 = inteval->stack[((hsi*1+7)*1+lsi)*1] + fp75;
LIBINT2_REALTYPE fp76;
fp76 = fp74 + fp77;
inteval->stack[((hsi*18+22)*1+lsi)*1] = fp76;
LIBINT2_REALTYPE fp71;
fp71 = inteval->AB_z[vi] * fp30;
LIBINT2_REALTYPE fp69;
fp69 = inteval->AB_y[vi] * inteval->stack[((hsi*1+6)*1+lsi)*1];
LIBINT2_REALTYPE fp68;
fp68 = inteval->stack[((hsi*1+8)*1+lsi)*1] + fp69;
LIBINT2_REALTYPE fp70;
fp70 = fp68 + fp71;
inteval->stack[((hsi*18+23)*1+lsi)*1] = fp70;
LIBINT2_REALTYPE fp27;
fp27 = inteval->AB_z[vi] * inteval->stack[((hsi*1+2)*1+lsi)*1];
LIBINT2_REALTYPE fp26;
fp26 = inteval->stack[((hsi*1+6)*1+lsi)*1] + fp27;
LIBINT2_REALTYPE fp67;
fp67 = inteval->AB_z[vi] * fp26;
LIBINT2_REALTYPE fp65;
fp65 = inteval->AB_z[vi] * inteval->stack[((hsi*1+6)*1+lsi)*1];
LIBINT2_REALTYPE fp64;
fp64 = inteval->stack[((hsi*1+9)*1+lsi)*1] + fp65;
LIBINT2_REALTYPE fp66;
fp66 = fp64 + fp67;
inteval->stack[((hsi*18+24)*1+lsi)*1] = fp66;
LIBINT2_REALTYPE fp25;
fp25 = inteval->AB_x[vi] * inteval->stack[((hsi*1+10)*1+lsi)*1];
LIBINT2_REALTYPE fp24;
fp24 = inteval->stack[((hsi*1+4)*1+lsi)*1] + fp25;
LIBINT2_REALTYPE fp55;
fp55 = inteval->AB_x[vi] * fp24;
LIBINT2_REALTYPE fp54;
fp54 = fp32 + fp55;
inteval->stack[((hsi*18+25)*1+lsi)*1] = fp54;
LIBINT2_REALTYPE fp53;
fp53 = inteval->AB_y[vi] * fp24;
LIBINT2_REALTYPE fp23;
fp23 = inteval->AB_x[vi] * inteval->stack[((hsi*1+11)*1+lsi)*1];
LIBINT2_REALTYPE fp22;
fp22 = inteval->stack[((hsi*1+7)*1+lsi)*1] + fp23;
LIBINT2_REALTYPE fp52;
fp52 = fp22 + fp53;
inteval->stack[((hsi*18+26)*1+lsi)*1] = fp52;
LIBINT2_REALTYPE fp19;
fp19 = inteval->AB_x[vi] * inteval->stack[((hsi*1+12)*1+lsi)*1];
LIBINT2_REALTYPE fp18;
fp18 = inteval->stack[((hsi*1+8)*1+lsi)*1] + fp19;
LIBINT2_REALTYPE fp73;
fp73 = inteval->AB_z[vi] * fp24;
LIBINT2_REALTYPE fp72;
fp72 = fp18 + fp73;
inteval->stack[((hsi*18+27)*1+lsi)*1] = fp72;
LIBINT2_REALTYPE fp57;
fp57 = inteval->AB_x[vi] * fp34;
LIBINT2_REALTYPE fp37;
fp37 = inteval->AB_x[vi] * inteval->stack[((hsi*1+1)*1+lsi)*1];
LIBINT2_REALTYPE fp36;
fp36 = inteval->stack[((hsi*1+0)*1+lsi)*1] + fp37;
LIBINT2_REALTYPE fp56;
fp56 = fp36 + fp57;
inteval->stack[((hsi*18+19)*1+lsi)*1] = fp56;
LIBINT2_REALTYPE fp13;
fp13 = inteval->AB_y[vi] * inteval->stack[((hsi*1+12)*1+lsi)*1];
LIBINT2_REALTYPE fp12;
fp12 = inteval->stack[((hsi*1+14)*1+lsi)*1] + fp13;
LIBINT2_REALTYPE fp21;
fp21 = inteval->AB_y[vi] * inteval->stack[((hsi*1+10)*1+lsi)*1];
LIBINT2_REALTYPE fp20;
fp20 = inteval->stack[((hsi*1+11)*1+lsi)*1] + fp21;
LIBINT2_REALTYPE fp63;
fp63 = inteval->AB_z[vi] * fp20;
LIBINT2_REALTYPE fp62;
fp62 = fp12 + fp63;
inteval->stack[((hsi*18+29)*1+lsi)*1] = fp62;
LIBINT2_REALTYPE fp17;
fp17 = inteval->AB_z[vi] * inteval->stack[((hsi*1+10)*1+lsi)*1];
LIBINT2_REALTYPE fp16;
fp16 = inteval->stack[((hsi*1+12)*1+lsi)*1] + fp17;
LIBINT2_REALTYPE fp61;
fp61 = inteval->AB_z[vi] * fp16;
LIBINT2_REALTYPE fp59;
fp59 = inteval->AB_z[vi] * inteval->stack[((hsi*1+12)*1+lsi)*1];
LIBINT2_REALTYPE fp58;
fp58 = inteval->stack[((hsi*1+15)*1+lsi)*1] + fp59;
LIBINT2_REALTYPE fp60;
fp60 = fp58 + fp61;
inteval->stack[((hsi*18+30)*1+lsi)*1] = fp60;
LIBINT2_REALTYPE fp11;
fp11 = inteval->AB_x[vi] * inteval->stack[((hsi*1+16)*1+lsi)*1];
LIBINT2_REALTYPE fp10;
fp10 = inteval->stack[((hsi*1+6)*1+lsi)*1] + fp11;
LIBINT2_REALTYPE fp49;
fp49 = inteval->AB_x[vi] * fp10;
LIBINT2_REALTYPE fp48;
fp48 = fp28 + fp49;
inteval->stack[((hsi*18+31)*1+lsi)*1] = fp48;
LIBINT2_REALTYPE fp47;
fp47 = inteval->AB_y[vi] * fp10;
LIBINT2_REALTYPE fp46;
fp46 = fp18 + fp47;
inteval->stack[((hsi*18+32)*1+lsi)*1] = fp46;
LIBINT2_REALTYPE fp45;
fp45 = inteval->AB_z[vi] * fp10;
LIBINT2_REALTYPE fp7;
fp7 = inteval->AB_x[vi] * inteval->stack[((hsi*1+17)*1+lsi)*1];
LIBINT2_REALTYPE fp6;
fp6 = inteval->stack[((hsi*1+9)*1+lsi)*1] + fp7;
LIBINT2_REALTYPE fp44;
fp44 = fp6 + fp45;
inteval->stack[((hsi*18+33)*1+lsi)*1] = fp44;
LIBINT2_REALTYPE fp9;
fp9 = inteval->AB_y[vi] * inteval->stack[((hsi*1+16)*1+lsi)*1];
LIBINT2_REALTYPE fp8;
fp8 = inteval->stack[((hsi*1+12)*1+lsi)*1] + fp9;
LIBINT2_REALTYPE fp43;
fp43 = inteval->AB_y[vi] * fp8;
LIBINT2_REALTYPE fp42;
fp42 = fp12 + fp43;
inteval->stack[((hsi*18+34)*1+lsi)*1] = fp42;
LIBINT2_REALTYPE fp41;
fp41 = inteval->AB_z[vi] * fp8;
LIBINT2_REALTYPE fp3;
fp3 = inteval->AB_y[vi] * inteval->stack[((hsi*1+17)*1+lsi)*1];
LIBINT2_REALTYPE fp2;
fp2 = inteval->stack[((hsi*1+15)*1+lsi)*1] + fp3;
LIBINT2_REALTYPE fp40;
fp40 = fp2 + fp41;
inteval->stack[((hsi*18+35)*1+lsi)*1] = fp40;
LIBINT2_REALTYPE fp5;
fp5 = inteval->AB_z[vi] * inteval->stack[((hsi*1+16)*1+lsi)*1];
LIBINT2_REALTYPE fp4;
fp4 = inteval->stack[((hsi*1+17)*1+lsi)*1] + fp5;
LIBINT2_REALTYPE fp39;
fp39 = inteval->AB_z[vi] * fp4;
LIBINT2_REALTYPE fp1;
fp1 = inteval->AB_z[vi] * inteval->stack[((hsi*1+17)*1+lsi)*1];
LIBINT2_REALTYPE fp0;
fp0 = inteval->stack[((hsi*1+18)*1+lsi)*1] + fp1;
LIBINT2_REALTYPE fp38;
fp38 = fp0 + fp39;
inteval->stack[((hsi*18+36)*1+lsi)*1] = fp38;
LIBINT2_REALTYPE fp51;
fp51 = inteval->AB_y[vi] * fp20;
LIBINT2_REALTYPE fp15;
fp15 = inteval->AB_y[vi] * inteval->stack[((hsi*1+11)*1+lsi)*1];
LIBINT2_REALTYPE fp14;
fp14 = inteval->stack[((hsi*1+13)*1+lsi)*1] + fp15;
LIBINT2_REALTYPE fp50;
fp50 = fp14 + fp51;
inteval->stack[((hsi*18+28)*1+lsi)*1] = fp50;
}
}
}
const int hsi = 0;
const int lsi = 0;
const int vi = 0;
inteval->targets[0] = &(inteval->stack[((hsi*18+19)*1+lsi)*1]);
/** Number of flops = 82 */
}

#ifdef __cplusplus
};
#endif