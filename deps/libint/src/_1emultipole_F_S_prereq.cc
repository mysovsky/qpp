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
void _1emultipole_F_S_prereq(const Libint_t* inteval, LIBINT2_REALTYPE* parent_stack) {

LIBINT2_REALTYPE*const  stack = parent_stack;
{
const int hsi = 0;
{
const int lsi = 0;
{
const int vi = 0;
LIBINT2_REALTYPE fp81;
fp81 = 0 + inteval->_0_Overlap_0_y[vi];
LIBINT2_REALTYPE fp135;
fp135 = 0 + fp81;
LIBINT2_REALTYPE fp95;
fp95 = 0 + inteval->_0_Overlap_0_x[vi];
LIBINT2_REALTYPE fp148;
fp148 = 0 + fp95;
LIBINT2_REALTYPE fp203;
fp203 = fp148 * fp135;
LIBINT2_REALTYPE fp102;
fp102 = 2 * inteval->oo2z[vi];
LIBINT2_REALTYPE fp60;
fp60 = 0 + inteval->_0_Overlap_0_z[vi];
LIBINT2_REALTYPE fp59;
fp59 = inteval->PA_z[vi] * fp60;
LIBINT2_REALTYPE fp54;
fp54 = fp102 * fp59;
LIBINT2_REALTYPE fp57;
fp57 = inteval->oo2z[vi] * fp60;
LIBINT2_REALTYPE fp58;
fp58 = inteval->PA_z[vi] * fp59;
LIBINT2_REALTYPE fp56;
fp56 = fp58 + fp57;
LIBINT2_REALTYPE fp55;
fp55 = inteval->PA_z[vi] * fp56;
LIBINT2_REALTYPE fp53;
fp53 = fp55 + fp54;
LIBINT2_REALTYPE fp116;
fp116 = 0 + fp53;
LIBINT2_REALTYPE fp106;
fp106 = inteval->BO_z[vi] * fp116;
LIBINT2_REALTYPE fp41;
fp41 = inteval->oo2z[vi] * fp56;
LIBINT2_REALTYPE fp52;
fp52 = inteval->PB_z[vi] * fp60;
LIBINT2_REALTYPE fp51;
fp51 = inteval->PA_z[vi] * fp52;
LIBINT2_REALTYPE fp50;
fp50 = fp51 + fp57;
LIBINT2_REALTYPE fp43;
fp43 = fp102 * fp50;
LIBINT2_REALTYPE fp46;
fp46 = inteval->oo2z[vi] * fp59;
LIBINT2_REALTYPE fp48;
fp48 = inteval->oo2z[vi] * fp52;
LIBINT2_REALTYPE fp49;
fp49 = inteval->PA_z[vi] * fp50;
LIBINT2_REALTYPE fp47;
fp47 = fp49 + fp48;
LIBINT2_REALTYPE fp45;
fp45 = fp47 + fp46;
LIBINT2_REALTYPE fp44;
fp44 = inteval->PA_z[vi] * fp45;
LIBINT2_REALTYPE fp42;
fp42 = fp44 + fp43;
LIBINT2_REALTYPE fp40;
fp40 = fp42 + fp41;
LIBINT2_REALTYPE fp104;
fp104 = 0 + fp40;
LIBINT2_REALTYPE fp105;
fp105 = fp104 + fp106;
LIBINT2_REALTYPE fp152;
fp152 = fp203 * fp105;
LIBINT2_REALTYPE fp39;
fp39 = fp152;
_libint2_static_api_inc1_short_(&(stack[((hsi*10+39)*1+lsi)*1]),&(fp39),1);
LIBINT2_REALTYPE fp117;
fp117 = 0 + fp56;
LIBINT2_REALTYPE fp109;
fp109 = inteval->BO_z[vi] * fp117;
LIBINT2_REALTYPE fp107;
fp107 = 0 + fp45;
LIBINT2_REALTYPE fp108;
fp108 = fp107 + fp109;
LIBINT2_REALTYPE fp80;
fp80 = inteval->PA_y[vi] * fp81;
LIBINT2_REALTYPE fp134;
fp134 = 0 + fp80;
LIBINT2_REALTYPE fp205;
fp205 = fp148 * fp134;
LIBINT2_REALTYPE fp153;
fp153 = fp205 * fp108;
LIBINT2_REALTYPE fp38;
fp38 = fp153;
_libint2_static_api_inc1_short_(&(stack[((hsi*10+38)*1+lsi)*1]),&(fp38),1);
LIBINT2_REALTYPE fp118;
fp118 = 0 + fp59;
LIBINT2_REALTYPE fp112;
fp112 = inteval->BO_z[vi] * fp118;
LIBINT2_REALTYPE fp110;
fp110 = 0 + fp50;
LIBINT2_REALTYPE fp111;
fp111 = fp110 + fp112;
LIBINT2_REALTYPE fp78;
fp78 = inteval->oo2z[vi] * fp81;
LIBINT2_REALTYPE fp79;
fp79 = inteval->PA_y[vi] * fp80;
LIBINT2_REALTYPE fp77;
fp77 = fp79 + fp78;
LIBINT2_REALTYPE fp133;
fp133 = 0 + fp77;
LIBINT2_REALTYPE fp207;
fp207 = fp148 * fp133;
LIBINT2_REALTYPE fp154;
fp154 = fp207 * fp111;
LIBINT2_REALTYPE fp37;
fp37 = fp154;
_libint2_static_api_inc1_short_(&(stack[((hsi*10+37)*1+lsi)*1]),&(fp37),1);
LIBINT2_REALTYPE fp119;
fp119 = 0 + fp60;
LIBINT2_REALTYPE fp115;
fp115 = inteval->BO_z[vi] * fp119;
LIBINT2_REALTYPE fp113;
fp113 = 0 + fp52;
LIBINT2_REALTYPE fp114;
fp114 = fp113 + fp115;
LIBINT2_REALTYPE fp75;
fp75 = fp102 * fp80;
LIBINT2_REALTYPE fp76;
fp76 = inteval->PA_y[vi] * fp77;
LIBINT2_REALTYPE fp74;
fp74 = fp76 + fp75;
LIBINT2_REALTYPE fp132;
fp132 = 0 + fp74;
LIBINT2_REALTYPE fp209;
fp209 = fp148 * fp132;
LIBINT2_REALTYPE fp155;
fp155 = fp209 * fp114;
LIBINT2_REALTYPE fp36;
fp36 = fp155;
_libint2_static_api_inc1_short_(&(stack[((hsi*10+36)*1+lsi)*1]),&(fp36),1);
LIBINT2_REALTYPE fp96;
fp96 = inteval->PA_x[vi] * fp95;
LIBINT2_REALTYPE fp149;
fp149 = 0 + fp96;
LIBINT2_REALTYPE fp211;
fp211 = fp149 * fp135;
LIBINT2_REALTYPE fp156;
fp156 = fp211 * fp108;
LIBINT2_REALTYPE fp35;
fp35 = fp156;
_libint2_static_api_inc1_short_(&(stack[((hsi*10+35)*1+lsi)*1]),&(fp35),1);
LIBINT2_REALTYPE fp213;
fp213 = fp149 * fp134;
LIBINT2_REALTYPE fp157;
fp157 = fp213 * fp111;
LIBINT2_REALTYPE fp34;
fp34 = fp157;
_libint2_static_api_inc1_short_(&(stack[((hsi*10+34)*1+lsi)*1]),&(fp34),1);
LIBINT2_REALTYPE fp215;
fp215 = fp149 * fp133;
LIBINT2_REALTYPE fp158;
fp158 = fp215 * fp114;
LIBINT2_REALTYPE fp33;
fp33 = fp158;
_libint2_static_api_inc1_short_(&(stack[((hsi*10+33)*1+lsi)*1]),&(fp33),1);
LIBINT2_REALTYPE fp98;
fp98 = inteval->oo2z[vi] * fp95;
LIBINT2_REALTYPE fp99;
fp99 = inteval->PA_x[vi] * fp96;
LIBINT2_REALTYPE fp97;
fp97 = fp99 + fp98;
LIBINT2_REALTYPE fp150;
fp150 = 0 + fp97;
LIBINT2_REALTYPE fp217;
fp217 = fp150 * fp135;
LIBINT2_REALTYPE fp159;
fp159 = fp217 * fp111;
LIBINT2_REALTYPE fp32;
fp32 = fp159;
_libint2_static_api_inc1_short_(&(stack[((hsi*10+32)*1+lsi)*1]),&(fp32),1);
LIBINT2_REALTYPE fp219;
fp219 = fp150 * fp134;
LIBINT2_REALTYPE fp160;
fp160 = fp219 * fp114;
LIBINT2_REALTYPE fp31;
fp31 = fp160;
_libint2_static_api_inc1_short_(&(stack[((hsi*10+31)*1+lsi)*1]),&(fp31),1);
LIBINT2_REALTYPE fp101;
fp101 = fp102 * fp96;
LIBINT2_REALTYPE fp103;
fp103 = inteval->PA_x[vi] * fp97;
LIBINT2_REALTYPE fp100;
fp100 = fp103 + fp101;
LIBINT2_REALTYPE fp151;
fp151 = 0 + fp100;
LIBINT2_REALTYPE fp221;
fp221 = fp151 * fp135;
LIBINT2_REALTYPE fp161;
fp161 = fp221 * fp114;
LIBINT2_REALTYPE fp30;
fp30 = fp161;
_libint2_static_api_inc1_short_(&(stack[((hsi*10+30)*1+lsi)*1]),&(fp30),1);
LIBINT2_REALTYPE fp131;
fp131 = inteval->BO_y[vi] * fp135;
LIBINT2_REALTYPE fp73;
fp73 = inteval->PB_y[vi] * fp81;
LIBINT2_REALTYPE fp129;
fp129 = 0 + fp73;
LIBINT2_REALTYPE fp130;
fp130 = fp129 + fp131;
LIBINT2_REALTYPE fp163;
fp163 = fp148 * fp130;
LIBINT2_REALTYPE fp162;
fp162 = fp163 * fp116;
LIBINT2_REALTYPE fp29;
fp29 = fp162;
_libint2_static_api_inc1_short_(&(stack[((hsi*10+29)*1+lsi)*1]),&(fp29),1);
LIBINT2_REALTYPE fp128;
fp128 = inteval->BO_y[vi] * fp134;
LIBINT2_REALTYPE fp72;
fp72 = inteval->PA_y[vi] * fp73;
LIBINT2_REALTYPE fp71;
fp71 = fp72 + fp78;
LIBINT2_REALTYPE fp126;
fp126 = 0 + fp71;
LIBINT2_REALTYPE fp127;
fp127 = fp126 + fp128;
LIBINT2_REALTYPE fp165;
fp165 = fp148 * fp127;
LIBINT2_REALTYPE fp164;
fp164 = fp165 * fp117;
LIBINT2_REALTYPE fp28;
fp28 = fp164;
_libint2_static_api_inc1_short_(&(stack[((hsi*10+28)*1+lsi)*1]),&(fp28),1);
LIBINT2_REALTYPE fp125;
fp125 = inteval->BO_y[vi] * fp133;
LIBINT2_REALTYPE fp67;
fp67 = inteval->oo2z[vi] * fp80;
LIBINT2_REALTYPE fp69;
fp69 = inteval->oo2z[vi] * fp73;
LIBINT2_REALTYPE fp70;
fp70 = inteval->PA_y[vi] * fp71;
LIBINT2_REALTYPE fp68;
fp68 = fp70 + fp69;
LIBINT2_REALTYPE fp66;
fp66 = fp68 + fp67;
LIBINT2_REALTYPE fp123;
fp123 = 0 + fp66;
LIBINT2_REALTYPE fp124;
fp124 = fp123 + fp125;
LIBINT2_REALTYPE fp167;
fp167 = fp148 * fp124;
LIBINT2_REALTYPE fp166;
fp166 = fp167 * fp118;
LIBINT2_REALTYPE fp27;
fp27 = fp166;
_libint2_static_api_inc1_short_(&(stack[((hsi*10+27)*1+lsi)*1]),&(fp27),1);
LIBINT2_REALTYPE fp122;
fp122 = inteval->BO_y[vi] * fp132;
LIBINT2_REALTYPE fp62;
fp62 = inteval->oo2z[vi] * fp77;
LIBINT2_REALTYPE fp64;
fp64 = fp102 * fp71;
LIBINT2_REALTYPE fp65;
fp65 = inteval->PA_y[vi] * fp66;
LIBINT2_REALTYPE fp63;
fp63 = fp65 + fp64;
LIBINT2_REALTYPE fp61;
fp61 = fp63 + fp62;
LIBINT2_REALTYPE fp120;
fp120 = 0 + fp61;
LIBINT2_REALTYPE fp121;
fp121 = fp120 + fp122;
LIBINT2_REALTYPE fp169;
fp169 = fp148 * fp121;
LIBINT2_REALTYPE fp168;
fp168 = fp169 * fp119;
LIBINT2_REALTYPE fp26;
fp26 = fp168;
_libint2_static_api_inc1_short_(&(stack[((hsi*10+26)*1+lsi)*1]),&(fp26),1);
LIBINT2_REALTYPE fp171;
fp171 = fp149 * fp130;
LIBINT2_REALTYPE fp170;
fp170 = fp171 * fp117;
LIBINT2_REALTYPE fp25;
fp25 = fp170;
_libint2_static_api_inc1_short_(&(stack[((hsi*10+25)*1+lsi)*1]),&(fp25),1);
LIBINT2_REALTYPE fp173;
fp173 = fp149 * fp127;
LIBINT2_REALTYPE fp172;
fp172 = fp173 * fp118;
LIBINT2_REALTYPE fp24;
fp24 = fp172;
_libint2_static_api_inc1_short_(&(stack[((hsi*10+24)*1+lsi)*1]),&(fp24),1);
LIBINT2_REALTYPE fp175;
fp175 = fp149 * fp124;
LIBINT2_REALTYPE fp174;
fp174 = fp175 * fp119;
LIBINT2_REALTYPE fp23;
fp23 = fp174;
_libint2_static_api_inc1_short_(&(stack[((hsi*10+23)*1+lsi)*1]),&(fp23),1);
LIBINT2_REALTYPE fp177;
fp177 = fp150 * fp130;
LIBINT2_REALTYPE fp176;
fp176 = fp177 * fp118;
LIBINT2_REALTYPE fp22;
fp22 = fp176;
_libint2_static_api_inc1_short_(&(stack[((hsi*10+22)*1+lsi)*1]),&(fp22),1);
LIBINT2_REALTYPE fp179;
fp179 = fp150 * fp127;
LIBINT2_REALTYPE fp178;
fp178 = fp179 * fp119;
LIBINT2_REALTYPE fp21;
fp21 = fp178;
_libint2_static_api_inc1_short_(&(stack[((hsi*10+21)*1+lsi)*1]),&(fp21),1);
LIBINT2_REALTYPE fp181;
fp181 = fp151 * fp130;
LIBINT2_REALTYPE fp180;
fp180 = fp181 * fp119;
LIBINT2_REALTYPE fp20;
fp20 = fp180;
_libint2_static_api_inc1_short_(&(stack[((hsi*10+20)*1+lsi)*1]),&(fp20),1);
LIBINT2_REALTYPE fp138;
fp138 = inteval->BO_x[vi] * fp148;
LIBINT2_REALTYPE fp82;
fp82 = inteval->PB_x[vi] * fp95;
LIBINT2_REALTYPE fp136;
fp136 = 0 + fp82;
LIBINT2_REALTYPE fp137;
fp137 = fp136 + fp138;
LIBINT2_REALTYPE fp183;
fp183 = fp137 * fp135;
LIBINT2_REALTYPE fp182;
fp182 = fp183 * fp116;
LIBINT2_REALTYPE fp19;
fp19 = fp182;
_libint2_static_api_inc1_short_(&(stack[((hsi*10+19)*1+lsi)*1]),&(fp19),1);
LIBINT2_REALTYPE fp185;
fp185 = fp137 * fp134;
LIBINT2_REALTYPE fp184;
fp184 = fp185 * fp117;
LIBINT2_REALTYPE fp18;
fp18 = fp184;
_libint2_static_api_inc1_short_(&(stack[((hsi*10+18)*1+lsi)*1]),&(fp18),1);
LIBINT2_REALTYPE fp187;
fp187 = fp137 * fp133;
LIBINT2_REALTYPE fp186;
fp186 = fp187 * fp118;
LIBINT2_REALTYPE fp17;
fp17 = fp186;
_libint2_static_api_inc1_short_(&(stack[((hsi*10+17)*1+lsi)*1]),&(fp17),1);
LIBINT2_REALTYPE fp189;
fp189 = fp137 * fp132;
LIBINT2_REALTYPE fp188;
fp188 = fp189 * fp119;
LIBINT2_REALTYPE fp16;
fp16 = fp188;
_libint2_static_api_inc1_short_(&(stack[((hsi*10+16)*1+lsi)*1]),&(fp16),1);
LIBINT2_REALTYPE fp141;
fp141 = inteval->BO_x[vi] * fp149;
LIBINT2_REALTYPE fp84;
fp84 = inteval->PA_x[vi] * fp82;
LIBINT2_REALTYPE fp83;
fp83 = fp84 + fp98;
LIBINT2_REALTYPE fp139;
fp139 = 0 + fp83;
LIBINT2_REALTYPE fp140;
fp140 = fp139 + fp141;
LIBINT2_REALTYPE fp191;
fp191 = fp140 * fp135;
LIBINT2_REALTYPE fp190;
fp190 = fp191 * fp117;
LIBINT2_REALTYPE fp15;
fp15 = fp190;
_libint2_static_api_inc1_short_(&(stack[((hsi*10+15)*1+lsi)*1]),&(fp15),1);
LIBINT2_REALTYPE fp193;
fp193 = fp140 * fp134;
LIBINT2_REALTYPE fp192;
fp192 = fp193 * fp118;
LIBINT2_REALTYPE fp14;
fp14 = fp192;
_libint2_static_api_inc1_short_(&(stack[((hsi*10+14)*1+lsi)*1]),&(fp14),1);
LIBINT2_REALTYPE fp195;
fp195 = fp140 * fp133;
LIBINT2_REALTYPE fp194;
fp194 = fp195 * fp119;
LIBINT2_REALTYPE fp13;
fp13 = fp194;
_libint2_static_api_inc1_short_(&(stack[((hsi*10+13)*1+lsi)*1]),&(fp13),1);
LIBINT2_REALTYPE fp144;
fp144 = inteval->BO_x[vi] * fp150;
LIBINT2_REALTYPE fp86;
fp86 = inteval->oo2z[vi] * fp96;
LIBINT2_REALTYPE fp88;
fp88 = inteval->oo2z[vi] * fp82;
LIBINT2_REALTYPE fp89;
fp89 = inteval->PA_x[vi] * fp83;
LIBINT2_REALTYPE fp87;
fp87 = fp89 + fp88;
LIBINT2_REALTYPE fp85;
fp85 = fp87 + fp86;
LIBINT2_REALTYPE fp142;
fp142 = 0 + fp85;
LIBINT2_REALTYPE fp143;
fp143 = fp142 + fp144;
LIBINT2_REALTYPE fp197;
fp197 = fp143 * fp135;
LIBINT2_REALTYPE fp196;
fp196 = fp197 * fp118;
LIBINT2_REALTYPE fp12;
fp12 = fp196;
_libint2_static_api_inc1_short_(&(stack[((hsi*10+12)*1+lsi)*1]),&(fp12),1);
LIBINT2_REALTYPE fp199;
fp199 = fp143 * fp134;
LIBINT2_REALTYPE fp198;
fp198 = fp199 * fp119;
LIBINT2_REALTYPE fp11;
fp11 = fp198;
_libint2_static_api_inc1_short_(&(stack[((hsi*10+11)*1+lsi)*1]),&(fp11),1);
LIBINT2_REALTYPE fp147;
fp147 = inteval->BO_x[vi] * fp151;
LIBINT2_REALTYPE fp91;
fp91 = inteval->oo2z[vi] * fp97;
LIBINT2_REALTYPE fp93;
fp93 = fp102 * fp83;
LIBINT2_REALTYPE fp94;
fp94 = inteval->PA_x[vi] * fp85;
LIBINT2_REALTYPE fp92;
fp92 = fp94 + fp93;
LIBINT2_REALTYPE fp90;
fp90 = fp92 + fp91;
LIBINT2_REALTYPE fp145;
fp145 = 0 + fp90;
LIBINT2_REALTYPE fp146;
fp146 = fp145 + fp147;
LIBINT2_REALTYPE fp201;
fp201 = fp146 * fp135;
LIBINT2_REALTYPE fp200;
fp200 = fp201 * fp119;
LIBINT2_REALTYPE fp10;
fp10 = fp200;
_libint2_static_api_inc1_short_(&(stack[((hsi*10+10)*1+lsi)*1]),&(fp10),1);
LIBINT2_REALTYPE fp202;
fp202 = fp203 * fp116;
LIBINT2_REALTYPE fp9;
fp9 = fp202;
_libint2_static_api_inc1_short_(&(stack[((hsi*10+9)*1+lsi)*1]),&(fp9),1);
LIBINT2_REALTYPE fp204;
fp204 = fp205 * fp117;
LIBINT2_REALTYPE fp8;
fp8 = fp204;
_libint2_static_api_inc1_short_(&(stack[((hsi*10+8)*1+lsi)*1]),&(fp8),1);
LIBINT2_REALTYPE fp206;
fp206 = fp207 * fp118;
LIBINT2_REALTYPE fp7;
fp7 = fp206;
_libint2_static_api_inc1_short_(&(stack[((hsi*10+7)*1+lsi)*1]),&(fp7),1);
LIBINT2_REALTYPE fp208;
fp208 = fp209 * fp119;
LIBINT2_REALTYPE fp6;
fp6 = fp208;
_libint2_static_api_inc1_short_(&(stack[((hsi*10+6)*1+lsi)*1]),&(fp6),1);
LIBINT2_REALTYPE fp210;
fp210 = fp211 * fp117;
LIBINT2_REALTYPE fp5;
fp5 = fp210;
_libint2_static_api_inc1_short_(&(stack[((hsi*10+5)*1+lsi)*1]),&(fp5),1);
LIBINT2_REALTYPE fp212;
fp212 = fp213 * fp118;
LIBINT2_REALTYPE fp4;
fp4 = fp212;
_libint2_static_api_inc1_short_(&(stack[((hsi*10+4)*1+lsi)*1]),&(fp4),1);
LIBINT2_REALTYPE fp214;
fp214 = fp215 * fp119;
LIBINT2_REALTYPE fp3;
fp3 = fp214;
_libint2_static_api_inc1_short_(&(stack[((hsi*10+3)*1+lsi)*1]),&(fp3),1);
LIBINT2_REALTYPE fp216;
fp216 = fp217 * fp118;
LIBINT2_REALTYPE fp2;
fp2 = fp216;
_libint2_static_api_inc1_short_(&(stack[((hsi*10+2)*1+lsi)*1]),&(fp2),1);
LIBINT2_REALTYPE fp218;
fp218 = fp219 * fp119;
LIBINT2_REALTYPE fp1;
fp1 = fp218;
_libint2_static_api_inc1_short_(&(stack[((hsi*10+1)*1+lsi)*1]),&(fp1),1);
LIBINT2_REALTYPE fp220;
fp220 = fp221 * fp119;
LIBINT2_REALTYPE fp0;
fp0 = fp220;
_libint2_static_api_inc1_short_(&(stack[((hsi*10+0)*1+lsi)*1]),&(fp0),1);
}
}
}
const int hsi = 0;
const int lsi = 0;
const int vi = 0;
/** Number of flops = 222 */
}

#ifdef __cplusplus
};
#endif
