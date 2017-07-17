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
#include <GenericContract.h>

#ifdef __cplusplus
extern "C" {
#endif
void _overlap_D_D_prereq(const Libint_t* inteval, LIBINT2_REALTYPE* parent_stack) {

LIBINT2_REALTYPE*const  stack = parent_stack;
{
const int hsi = 0;
{
const int lsi = 0;
{
const int vi = 0;
LIBINT2_REALTYPE fp102;
fp102 = 2 * inteval->oo2z[vi];
LIBINT2_REALTYPE fp58;
fp58 = 0 + inteval->_0_Overlap_0_z[vi];
LIBINT2_REALTYPE fp57;
fp57 = inteval->PB_z[vi] * fp58;
LIBINT2_REALTYPE fp49;
fp49 = fp102 * fp57;
LIBINT2_REALTYPE fp55;
fp55 = inteval->oo2z[vi] * fp58;
LIBINT2_REALTYPE fp56;
fp56 = inteval->PB_z[vi] * fp57;
LIBINT2_REALTYPE fp54;
fp54 = fp56 + fp55;
LIBINT2_REALTYPE fp50;
fp50 = inteval->PA_z[vi] * fp54;
LIBINT2_REALTYPE fp48;
fp48 = fp50 + fp49;
LIBINT2_REALTYPE fp109;
fp109 = 0 + fp48;
LIBINT2_REALTYPE fp81;
fp81 = 0 + inteval->_0_Overlap_0_y[vi];
LIBINT2_REALTYPE fp123;
fp123 = 0 + fp81;
LIBINT2_REALTYPE fp82;
fp82 = 0 + inteval->_0_Overlap_0_x[vi];
LIBINT2_REALTYPE fp83;
fp83 = inteval->PA_x[vi] * fp82;
LIBINT2_REALTYPE fp127;
fp127 = 0 + fp83;
LIBINT2_REALTYPE fp170;
fp170 = fp127 * fp123;
LIBINT2_REALTYPE fp169;
fp169 = fp170 * fp109;
LIBINT2_REALTYPE fp17;
fp17 = fp169;
_libint2_static_api_inc1_short_(&(stack[((hsi*36+17)*1+lsi)*1]),&(fp17),1);
LIBINT2_REALTYPE fp52;
fp52 = inteval->PA_z[vi] * fp57;
LIBINT2_REALTYPE fp51;
fp51 = fp52 + fp55;
LIBINT2_REALTYPE fp110;
fp110 = 0 + fp51;
LIBINT2_REALTYPE fp80;
fp80 = inteval->PB_y[vi] * fp81;
LIBINT2_REALTYPE fp122;
fp122 = 0 + fp80;
LIBINT2_REALTYPE fp172;
fp172 = fp127 * fp122;
LIBINT2_REALTYPE fp171;
fp171 = fp172 * fp110;
LIBINT2_REALTYPE fp16;
fp16 = fp171;
_libint2_static_api_inc1_short_(&(stack[((hsi*36+16)*1+lsi)*1]),&(fp16),1);
LIBINT2_REALTYPE fp53;
fp53 = inteval->PA_z[vi] * fp58;
LIBINT2_REALTYPE fp111;
fp111 = 0 + fp53;
LIBINT2_REALTYPE fp78;
fp78 = inteval->oo2z[vi] * fp81;
LIBINT2_REALTYPE fp79;
fp79 = inteval->PB_y[vi] * fp80;
LIBINT2_REALTYPE fp77;
fp77 = fp79 + fp78;
LIBINT2_REALTYPE fp121;
fp121 = 0 + fp77;
LIBINT2_REALTYPE fp174;
fp174 = fp127 * fp121;
LIBINT2_REALTYPE fp173;
fp173 = fp174 * fp111;
LIBINT2_REALTYPE fp15;
fp15 = fp173;
_libint2_static_api_inc1_short_(&(stack[((hsi*36+15)*1+lsi)*1]),&(fp15),1);
LIBINT2_REALTYPE fp95;
fp95 = inteval->oo2z[vi] * fp82;
LIBINT2_REALTYPE fp84;
fp84 = inteval->PB_x[vi] * fp82;
LIBINT2_REALTYPE fp93;
fp93 = inteval->PA_x[vi] * fp84;
LIBINT2_REALTYPE fp92;
fp92 = fp93 + fp95;
LIBINT2_REALTYPE fp128;
fp128 = 0 + fp92;
LIBINT2_REALTYPE fp176;
fp176 = fp128 * fp123;
LIBINT2_REALTYPE fp175;
fp175 = fp176 * fp110;
LIBINT2_REALTYPE fp14;
fp14 = fp175;
_libint2_static_api_inc1_short_(&(stack[((hsi*36+14)*1+lsi)*1]),&(fp14),1);
LIBINT2_REALTYPE fp178;
fp178 = fp128 * fp122;
LIBINT2_REALTYPE fp177;
fp177 = fp178 * fp111;
LIBINT2_REALTYPE fp13;
fp13 = fp177;
_libint2_static_api_inc1_short_(&(stack[((hsi*36+13)*1+lsi)*1]),&(fp13),1);
LIBINT2_REALTYPE fp98;
fp98 = fp102 * fp84;
LIBINT2_REALTYPE fp96;
fp96 = inteval->PB_x[vi] * fp84;
LIBINT2_REALTYPE fp94;
fp94 = fp96 + fp95;
LIBINT2_REALTYPE fp99;
fp99 = inteval->PA_x[vi] * fp94;
LIBINT2_REALTYPE fp97;
fp97 = fp99 + fp98;
LIBINT2_REALTYPE fp129;
fp129 = 0 + fp97;
LIBINT2_REALTYPE fp180;
fp180 = fp129 * fp123;
LIBINT2_REALTYPE fp179;
fp179 = fp180 * fp111;
LIBINT2_REALTYPE fp12;
fp12 = fp179;
_libint2_static_api_inc1_short_(&(stack[((hsi*36+12)*1+lsi)*1]),&(fp12),1);
LIBINT2_REALTYPE fp112;
fp112 = 0 + fp54;
LIBINT2_REALTYPE fp76;
fp76 = inteval->PA_y[vi] * fp81;
LIBINT2_REALTYPE fp120;
fp120 = 0 + fp76;
LIBINT2_REALTYPE fp182;
fp182 = fp127 * fp120;
LIBINT2_REALTYPE fp181;
fp181 = fp182 * fp112;
LIBINT2_REALTYPE fp11;
fp11 = fp181;
_libint2_static_api_inc1_short_(&(stack[((hsi*36+11)*1+lsi)*1]),&(fp11),1);
LIBINT2_REALTYPE fp113;
fp113 = 0 + fp57;
LIBINT2_REALTYPE fp75;
fp75 = inteval->PA_y[vi] * fp80;
LIBINT2_REALTYPE fp74;
fp74 = fp75 + fp78;
LIBINT2_REALTYPE fp119;
fp119 = 0 + fp74;
LIBINT2_REALTYPE fp184;
fp184 = fp127 * fp119;
LIBINT2_REALTYPE fp183;
fp183 = fp184 * fp113;
LIBINT2_REALTYPE fp10;
fp10 = fp183;
_libint2_static_api_inc1_short_(&(stack[((hsi*36+10)*1+lsi)*1]),&(fp10),1);
LIBINT2_REALTYPE fp114;
fp114 = 0 + fp58;
LIBINT2_REALTYPE fp101;
fp101 = fp102 * fp92;
LIBINT2_REALTYPE fp104;
fp104 = inteval->oo2z[vi] * fp94;
LIBINT2_REALTYPE fp105;
fp105 = inteval->PA_x[vi] * fp97;
LIBINT2_REALTYPE fp103;
fp103 = fp105 + fp104;
LIBINT2_REALTYPE fp100;
fp100 = fp103 + fp101;
LIBINT2_REALTYPE fp132;
fp132 = 0 + fp100;
LIBINT2_REALTYPE fp204;
fp204 = fp132 * fp123;
LIBINT2_REALTYPE fp203;
fp203 = fp204 * fp114;
LIBINT2_REALTYPE fp0;
fp0 = fp203;
_libint2_static_api_inc1_short_(&(stack[((hsi*36+0)*1+lsi)*1]),&(fp0),1);
LIBINT2_REALTYPE fp188;
fp188 = fp128 * fp120;
LIBINT2_REALTYPE fp187;
fp187 = fp188 * fp113;
LIBINT2_REALTYPE fp8;
fp8 = fp187;
_libint2_static_api_inc1_short_(&(stack[((hsi*36+8)*1+lsi)*1]),&(fp8),1);
LIBINT2_REALTYPE fp190;
fp190 = fp128 * fp119;
LIBINT2_REALTYPE fp189;
fp189 = fp190 * fp114;
LIBINT2_REALTYPE fp7;
fp7 = fp189;
_libint2_static_api_inc1_short_(&(stack[((hsi*36+7)*1+lsi)*1]),&(fp7),1);
LIBINT2_REALTYPE fp192;
fp192 = fp129 * fp120;
LIBINT2_REALTYPE fp191;
fp191 = fp192 * fp114;
LIBINT2_REALTYPE fp6;
fp6 = fp191;
_libint2_static_api_inc1_short_(&(stack[((hsi*36+6)*1+lsi)*1]),&(fp6),1);
LIBINT2_REALTYPE fp86;
fp86 = inteval->PA_x[vi] * fp83;
LIBINT2_REALTYPE fp85;
fp85 = fp86 + fp95;
LIBINT2_REALTYPE fp130;
fp130 = 0 + fp85;
LIBINT2_REALTYPE fp194;
fp194 = fp130 * fp123;
LIBINT2_REALTYPE fp193;
fp193 = fp194 * fp112;
LIBINT2_REALTYPE fp5;
fp5 = fp193;
_libint2_static_api_inc1_short_(&(stack[((hsi*36+5)*1+lsi)*1]),&(fp5),1);
LIBINT2_REALTYPE fp196;
fp196 = fp130 * fp122;
LIBINT2_REALTYPE fp195;
fp195 = fp196 * fp113;
LIBINT2_REALTYPE fp4;
fp4 = fp195;
_libint2_static_api_inc1_short_(&(stack[((hsi*36+4)*1+lsi)*1]),&(fp4),1);
LIBINT2_REALTYPE fp198;
fp198 = fp130 * fp121;
LIBINT2_REALTYPE fp197;
fp197 = fp198 * fp114;
LIBINT2_REALTYPE fp3;
fp3 = fp197;
_libint2_static_api_inc1_short_(&(stack[((hsi*36+3)*1+lsi)*1]),&(fp3),1);
LIBINT2_REALTYPE fp88;
fp88 = inteval->oo2z[vi] * fp83;
LIBINT2_REALTYPE fp90;
fp90 = inteval->oo2z[vi] * fp84;
LIBINT2_REALTYPE fp91;
fp91 = inteval->PA_x[vi] * fp92;
LIBINT2_REALTYPE fp89;
fp89 = fp91 + fp90;
LIBINT2_REALTYPE fp87;
fp87 = fp89 + fp88;
LIBINT2_REALTYPE fp131;
fp131 = 0 + fp87;
LIBINT2_REALTYPE fp200;
fp200 = fp131 * fp123;
LIBINT2_REALTYPE fp199;
fp199 = fp200 * fp113;
LIBINT2_REALTYPE fp2;
fp2 = fp199;
_libint2_static_api_inc1_short_(&(stack[((hsi*36+2)*1+lsi)*1]),&(fp2),1);
LIBINT2_REALTYPE fp202;
fp202 = fp131 * fp122;
LIBINT2_REALTYPE fp201;
fp201 = fp202 * fp114;
LIBINT2_REALTYPE fp1;
fp1 = fp201;
_libint2_static_api_inc1_short_(&(stack[((hsi*36+1)*1+lsi)*1]),&(fp1),1);
LIBINT2_REALTYPE fp72;
fp72 = fp102 * fp80;
LIBINT2_REALTYPE fp73;
fp73 = inteval->PA_y[vi] * fp77;
LIBINT2_REALTYPE fp71;
fp71 = fp73 + fp72;
LIBINT2_REALTYPE fp118;
fp118 = 0 + fp71;
LIBINT2_REALTYPE fp186;
fp186 = fp127 * fp118;
LIBINT2_REALTYPE fp185;
fp185 = fp186 * fp114;
LIBINT2_REALTYPE fp9;
fp9 = fp185;
_libint2_static_api_inc1_short_(&(stack[((hsi*36+9)*1+lsi)*1]),&(fp9),1);
LIBINT2_REALTYPE fp37;
fp37 = fp102 * fp51;
LIBINT2_REALTYPE fp39;
fp39 = inteval->oo2z[vi] * fp54;
LIBINT2_REALTYPE fp40;
fp40 = inteval->PA_z[vi] * fp48;
LIBINT2_REALTYPE fp38;
fp38 = fp40 + fp39;
LIBINT2_REALTYPE fp36;
fp36 = fp38 + fp37;
LIBINT2_REALTYPE fp106;
fp106 = 0 + fp36;
LIBINT2_REALTYPE fp124;
fp124 = 0 + fp82;
LIBINT2_REALTYPE fp134;
fp134 = fp124 * fp123;
LIBINT2_REALTYPE fp133;
fp133 = fp134 * fp106;
LIBINT2_REALTYPE fp35;
fp35 = fp133;
_libint2_static_api_inc1_short_(&(stack[((hsi*36+35)*1+lsi)*1]),&(fp35),1);
LIBINT2_REALTYPE fp42;
fp42 = inteval->oo2z[vi] * fp53;
LIBINT2_REALTYPE fp44;
fp44 = inteval->oo2z[vi] * fp57;
LIBINT2_REALTYPE fp45;
fp45 = inteval->PA_z[vi] * fp51;
LIBINT2_REALTYPE fp43;
fp43 = fp45 + fp44;
LIBINT2_REALTYPE fp41;
fp41 = fp43 + fp42;
LIBINT2_REALTYPE fp107;
fp107 = 0 + fp41;
LIBINT2_REALTYPE fp136;
fp136 = fp124 * fp122;
LIBINT2_REALTYPE fp135;
fp135 = fp136 * fp107;
LIBINT2_REALTYPE fp34;
fp34 = fp135;
_libint2_static_api_inc1_short_(&(stack[((hsi*36+34)*1+lsi)*1]),&(fp34),1);
LIBINT2_REALTYPE fp47;
fp47 = inteval->PA_z[vi] * fp53;
LIBINT2_REALTYPE fp46;
fp46 = fp47 + fp55;
LIBINT2_REALTYPE fp108;
fp108 = 0 + fp46;
LIBINT2_REALTYPE fp138;
fp138 = fp124 * fp121;
LIBINT2_REALTYPE fp137;
fp137 = fp138 * fp108;
LIBINT2_REALTYPE fp33;
fp33 = fp137;
_libint2_static_api_inc1_short_(&(stack[((hsi*36+33)*1+lsi)*1]),&(fp33),1);
LIBINT2_REALTYPE fp125;
fp125 = 0 + fp84;
LIBINT2_REALTYPE fp140;
fp140 = fp125 * fp123;
LIBINT2_REALTYPE fp139;
fp139 = fp140 * fp107;
LIBINT2_REALTYPE fp32;
fp32 = fp139;
_libint2_static_api_inc1_short_(&(stack[((hsi*36+32)*1+lsi)*1]),&(fp32),1);
LIBINT2_REALTYPE fp142;
fp142 = fp125 * fp122;
LIBINT2_REALTYPE fp141;
fp141 = fp142 * fp108;
LIBINT2_REALTYPE fp31;
fp31 = fp141;
_libint2_static_api_inc1_short_(&(stack[((hsi*36+31)*1+lsi)*1]),&(fp31),1);
LIBINT2_REALTYPE fp126;
fp126 = 0 + fp94;
LIBINT2_REALTYPE fp144;
fp144 = fp126 * fp123;
LIBINT2_REALTYPE fp143;
fp143 = fp144 * fp108;
LIBINT2_REALTYPE fp30;
fp30 = fp143;
_libint2_static_api_inc1_short_(&(stack[((hsi*36+30)*1+lsi)*1]),&(fp30),1);
LIBINT2_REALTYPE fp146;
fp146 = fp124 * fp120;
LIBINT2_REALTYPE fp145;
fp145 = fp146 * fp109;
LIBINT2_REALTYPE fp29;
fp29 = fp145;
_libint2_static_api_inc1_short_(&(stack[((hsi*36+29)*1+lsi)*1]),&(fp29),1);
LIBINT2_REALTYPE fp148;
fp148 = fp124 * fp119;
LIBINT2_REALTYPE fp147;
fp147 = fp148 * fp110;
LIBINT2_REALTYPE fp28;
fp28 = fp147;
_libint2_static_api_inc1_short_(&(stack[((hsi*36+28)*1+lsi)*1]),&(fp28),1);
LIBINT2_REALTYPE fp70;
fp70 = inteval->PA_y[vi] * fp76;
LIBINT2_REALTYPE fp69;
fp69 = fp70 + fp78;
LIBINT2_REALTYPE fp117;
fp117 = 0 + fp69;
LIBINT2_REALTYPE fp168;
fp168 = fp126 * fp117;
LIBINT2_REALTYPE fp167;
fp167 = fp168 * fp114;
LIBINT2_REALTYPE fp18;
fp18 = fp167;
_libint2_static_api_inc1_short_(&(stack[((hsi*36+18)*1+lsi)*1]),&(fp18),1);
LIBINT2_REALTYPE fp152;
fp152 = fp125 * fp120;
LIBINT2_REALTYPE fp151;
fp151 = fp152 * fp110;
LIBINT2_REALTYPE fp26;
fp26 = fp151;
_libint2_static_api_inc1_short_(&(stack[((hsi*36+26)*1+lsi)*1]),&(fp26),1);
LIBINT2_REALTYPE fp154;
fp154 = fp125 * fp119;
LIBINT2_REALTYPE fp153;
fp153 = fp154 * fp111;
LIBINT2_REALTYPE fp25;
fp25 = fp153;
_libint2_static_api_inc1_short_(&(stack[((hsi*36+25)*1+lsi)*1]),&(fp25),1);
LIBINT2_REALTYPE fp156;
fp156 = fp126 * fp120;
LIBINT2_REALTYPE fp155;
fp155 = fp156 * fp111;
LIBINT2_REALTYPE fp24;
fp24 = fp155;
_libint2_static_api_inc1_short_(&(stack[((hsi*36+24)*1+lsi)*1]),&(fp24),1);
LIBINT2_REALTYPE fp158;
fp158 = fp124 * fp117;
LIBINT2_REALTYPE fp157;
fp157 = fp158 * fp112;
LIBINT2_REALTYPE fp23;
fp23 = fp157;
_libint2_static_api_inc1_short_(&(stack[((hsi*36+23)*1+lsi)*1]),&(fp23),1);
LIBINT2_REALTYPE fp65;
fp65 = inteval->oo2z[vi] * fp76;
LIBINT2_REALTYPE fp67;
fp67 = inteval->oo2z[vi] * fp80;
LIBINT2_REALTYPE fp68;
fp68 = inteval->PA_y[vi] * fp74;
LIBINT2_REALTYPE fp66;
fp66 = fp68 + fp67;
LIBINT2_REALTYPE fp64;
fp64 = fp66 + fp65;
LIBINT2_REALTYPE fp116;
fp116 = 0 + fp64;
LIBINT2_REALTYPE fp160;
fp160 = fp124 * fp116;
LIBINT2_REALTYPE fp159;
fp159 = fp160 * fp113;
LIBINT2_REALTYPE fp22;
fp22 = fp159;
_libint2_static_api_inc1_short_(&(stack[((hsi*36+22)*1+lsi)*1]),&(fp22),1);
LIBINT2_REALTYPE fp60;
fp60 = fp102 * fp74;
LIBINT2_REALTYPE fp62;
fp62 = inteval->oo2z[vi] * fp77;
LIBINT2_REALTYPE fp63;
fp63 = inteval->PA_y[vi] * fp71;
LIBINT2_REALTYPE fp61;
fp61 = fp63 + fp62;
LIBINT2_REALTYPE fp59;
fp59 = fp61 + fp60;
LIBINT2_REALTYPE fp115;
fp115 = 0 + fp59;
LIBINT2_REALTYPE fp162;
fp162 = fp124 * fp115;
LIBINT2_REALTYPE fp161;
fp161 = fp162 * fp114;
LIBINT2_REALTYPE fp21;
fp21 = fp161;
_libint2_static_api_inc1_short_(&(stack[((hsi*36+21)*1+lsi)*1]),&(fp21),1);
LIBINT2_REALTYPE fp164;
fp164 = fp125 * fp117;
LIBINT2_REALTYPE fp163;
fp163 = fp164 * fp113;
LIBINT2_REALTYPE fp20;
fp20 = fp163;
_libint2_static_api_inc1_short_(&(stack[((hsi*36+20)*1+lsi)*1]),&(fp20),1);
LIBINT2_REALTYPE fp166;
fp166 = fp125 * fp116;
LIBINT2_REALTYPE fp165;
fp165 = fp166 * fp114;
LIBINT2_REALTYPE fp19;
fp19 = fp165;
_libint2_static_api_inc1_short_(&(stack[((hsi*36+19)*1+lsi)*1]),&(fp19),1);
LIBINT2_REALTYPE fp150;
fp150 = fp124 * fp118;
LIBINT2_REALTYPE fp149;
fp149 = fp150 * fp111;
LIBINT2_REALTYPE fp27;
fp27 = fp149;
_libint2_static_api_inc1_short_(&(stack[((hsi*36+27)*1+lsi)*1]),&(fp27),1);
}
}
}
const int hsi = 0;
const int lsi = 0;
const int vi = 0;
/** Number of flops = 205 */
}

#ifdef __cplusplus
};
#endif
