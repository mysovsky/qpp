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
void _aB_P__0__P__1___TwoPRep_S__0__P__1___Ab__up_0_prereq(const Libint_t* inteval, LIBINT2_REALTYPE* parent_stack) {

LIBINT2_REALTYPE*const  stack = parent_stack;
{
const int hsi = 0;
{
const int lsi = 0;
{
const int vi = 0;
LIBINT2_REALTYPE fp136;
fp136 = inteval->WQ_z[vi] * inteval->_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_2[vi];
LIBINT2_REALTYPE fp137;
fp137 = inteval->QC_z[vi] * inteval->_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_1[vi];
LIBINT2_REALTYPE fp135;
fp135 = fp137 + fp136;
LIBINT2_REALTYPE fp193;
fp193 = 2 * inteval->oo2ze[vi];
LIBINT2_REALTYPE fp28;
fp28 = fp193 * fp135;
LIBINT2_REALTYPE fp187;
fp187 = 1 * inteval->oo2e[vi];
LIBINT2_REALTYPE fp179;
fp179 = inteval->roe[vi] * inteval->_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_2[vi];
LIBINT2_REALTYPE fp178;
fp178 = inteval->_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_1[vi] - fp179;
LIBINT2_REALTYPE fp177;
fp177 = fp187 * fp178;
LIBINT2_REALTYPE fp107;
fp107 = inteval->WQ_z[vi] * inteval->_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_3[vi];
LIBINT2_REALTYPE fp108;
fp108 = inteval->QC_z[vi] * inteval->_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_2[vi];
LIBINT2_REALTYPE fp106;
fp106 = fp108 + fp107;
LIBINT2_REALTYPE fp97;
fp97 = inteval->WQ_z[vi] * fp106;
LIBINT2_REALTYPE fp98;
fp98 = inteval->QC_z[vi] * fp135;
LIBINT2_REALTYPE fp96;
fp96 = fp98 + fp97;
LIBINT2_REALTYPE fp95;
fp95 = fp96 + fp177;
LIBINT2_REALTYPE fp30;
fp30 = inteval->WP_z[vi] * fp95;
LIBINT2_REALTYPE fp186;
fp186 = inteval->roe[vi] * inteval->_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_1[vi];
LIBINT2_REALTYPE fp185;
fp185 = inteval->_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_0[vi] - fp186;
LIBINT2_REALTYPE fp184;
fp184 = fp187 * fp185;
LIBINT2_REALTYPE fp101;
fp101 = inteval->WQ_z[vi] * fp135;
LIBINT2_REALTYPE fp110;
fp110 = inteval->WQ_z[vi] * inteval->_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_1[vi];
LIBINT2_REALTYPE fp111;
fp111 = inteval->QC_z[vi] * inteval->_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_0[vi];
LIBINT2_REALTYPE fp109;
fp109 = fp111 + fp110;
LIBINT2_REALTYPE fp102;
fp102 = inteval->QC_z[vi] * fp109;
LIBINT2_REALTYPE fp100;
fp100 = fp102 + fp101;
LIBINT2_REALTYPE fp99;
fp99 = fp100 + fp184;
LIBINT2_REALTYPE fp31;
fp31 = inteval->PA_z[vi] * fp99;
LIBINT2_REALTYPE fp29;
fp29 = fp31 + fp30;
LIBINT2_REALTYPE fp27;
fp27 = fp29 + fp28;
LIBINT2_REALTYPE fp26;
fp26 = fp27;
_libint2_static_api_inc1_short_(&(stack[((hsi*1+26)*1+lsi)*1]),&(fp26),1);
LIBINT2_REALTYPE fp169;
fp169 = 1 * inteval->oo2ze[vi];
LIBINT2_REALTYPE fp168;
fp168 = fp169 * inteval->_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_1[vi];
LIBINT2_REALTYPE fp44;
fp44 = inteval->WP_z[vi] * fp135;
LIBINT2_REALTYPE fp45;
fp45 = inteval->PA_z[vi] * fp109;
LIBINT2_REALTYPE fp43;
fp43 = fp45 + fp44;
LIBINT2_REALTYPE fp42;
fp42 = fp43 + fp168;
LIBINT2_REALTYPE fp23;
fp23 = fp42;
_libint2_static_api_inc1_short_(&(stack[((hsi*1+25)*1+lsi)*1]),&(fp23),1);
LIBINT2_REALTYPE fp163;
fp163 = fp169 * inteval->_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_2[vi];
LIBINT2_REALTYPE fp40;
fp40 = inteval->WP_z[vi] * fp106;
LIBINT2_REALTYPE fp41;
fp41 = inteval->PA_z[vi] * fp135;
LIBINT2_REALTYPE fp39;
fp39 = fp41 + fp40;
LIBINT2_REALTYPE fp38;
fp38 = fp39 + fp163;
LIBINT2_REALTYPE fp33;
fp33 = inteval->WQ_y[vi] * fp38;
LIBINT2_REALTYPE fp34;
fp34 = inteval->QC_y[vi] * fp23;
LIBINT2_REALTYPE fp32;
fp32 = fp34 + fp33;
LIBINT2_REALTYPE fp25;
fp25 = fp32;
_libint2_static_api_inc1_short_(&(stack[((hsi*1+24)*1+lsi)*1]),&(fp25),1);
LIBINT2_REALTYPE fp116;
fp116 = inteval->WQ_y[vi] * inteval->_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_3[vi];
LIBINT2_REALTYPE fp117;
fp117 = inteval->QC_y[vi] * inteval->_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_2[vi];
LIBINT2_REALTYPE fp115;
fp115 = fp117 + fp116;
LIBINT2_REALTYPE fp120;
fp120 = inteval->WQ_y[vi] * fp115;
LIBINT2_REALTYPE fp145;
fp145 = inteval->WQ_y[vi] * inteval->_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_2[vi];
LIBINT2_REALTYPE fp146;
fp146 = inteval->QC_y[vi] * inteval->_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_1[vi];
LIBINT2_REALTYPE fp144;
fp144 = fp146 + fp145;
LIBINT2_REALTYPE fp121;
fp121 = inteval->QC_y[vi] * fp144;
LIBINT2_REALTYPE fp119;
fp119 = fp121 + fp120;
LIBINT2_REALTYPE fp118;
fp118 = fp119 + fp177;
LIBINT2_REALTYPE fp36;
fp36 = inteval->WP_z[vi] * fp118;
LIBINT2_REALTYPE fp124;
fp124 = inteval->WQ_y[vi] * fp144;
LIBINT2_REALTYPE fp130;
fp130 = inteval->WQ_y[vi] * inteval->_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_1[vi];
LIBINT2_REALTYPE fp131;
fp131 = inteval->QC_y[vi] * inteval->_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_0[vi];
LIBINT2_REALTYPE fp129;
fp129 = fp131 + fp130;
LIBINT2_REALTYPE fp125;
fp125 = inteval->QC_y[vi] * fp129;
LIBINT2_REALTYPE fp123;
fp123 = fp125 + fp124;
LIBINT2_REALTYPE fp122;
fp122 = fp123 + fp184;
LIBINT2_REALTYPE fp37;
fp37 = inteval->PA_z[vi] * fp122;
LIBINT2_REALTYPE fp35;
fp35 = fp37 + fp36;
LIBINT2_REALTYPE fp24;
fp24 = fp35;
_libint2_static_api_inc1_short_(&(stack[((hsi*1+23)*1+lsi)*1]),&(fp24),1);
LIBINT2_REALTYPE fp50;
fp50 = inteval->WP_z[vi] * fp144;
LIBINT2_REALTYPE fp51;
fp51 = inteval->PA_z[vi] * fp129;
LIBINT2_REALTYPE fp49;
fp49 = fp51 + fp50;
LIBINT2_REALTYPE fp21;
fp21 = fp49;
_libint2_static_api_inc1_short_(&(stack[((hsi*1+22)*1+lsi)*1]),&(fp21),1);
LIBINT2_REALTYPE fp47;
fp47 = inteval->WQ_x[vi] * fp38;
LIBINT2_REALTYPE fp48;
fp48 = inteval->QC_x[vi] * fp23;
LIBINT2_REALTYPE fp46;
fp46 = fp48 + fp47;
LIBINT2_REALTYPE fp22;
fp22 = fp46;
_libint2_static_api_inc1_short_(&(stack[((hsi*1+21)*1+lsi)*1]),&(fp22),1);
LIBINT2_REALTYPE fp157;
fp157 = inteval->WQ_x[vi] * inteval->_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_3[vi];
LIBINT2_REALTYPE fp158;
fp158 = inteval->QC_x[vi] * inteval->_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_2[vi];
LIBINT2_REALTYPE fp156;
fp156 = fp158 + fp157;
LIBINT2_REALTYPE fp148;
fp148 = inteval->WQ_y[vi] * fp156;
LIBINT2_REALTYPE fp174;
fp174 = inteval->WQ_x[vi] * inteval->_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_2[vi];
LIBINT2_REALTYPE fp175;
fp175 = inteval->QC_x[vi] * inteval->_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_1[vi];
LIBINT2_REALTYPE fp173;
fp173 = fp175 + fp174;
LIBINT2_REALTYPE fp149;
fp149 = inteval->QC_y[vi] * fp173;
LIBINT2_REALTYPE fp147;
fp147 = fp149 + fp148;
LIBINT2_REALTYPE fp53;
fp53 = inteval->WP_z[vi] * fp147;
LIBINT2_REALTYPE fp151;
fp151 = inteval->WQ_y[vi] * fp173;
LIBINT2_REALTYPE fp160;
fp160 = inteval->WQ_x[vi] * inteval->_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_1[vi];
LIBINT2_REALTYPE fp161;
fp161 = inteval->QC_x[vi] * inteval->_aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_0[vi];
LIBINT2_REALTYPE fp159;
fp159 = fp161 + fp160;
LIBINT2_REALTYPE fp152;
fp152 = inteval->QC_y[vi] * fp159;
LIBINT2_REALTYPE fp150;
fp150 = fp152 + fp151;
LIBINT2_REALTYPE fp54;
fp54 = inteval->PA_z[vi] * fp150;
LIBINT2_REALTYPE fp52;
fp52 = fp54 + fp53;
LIBINT2_REALTYPE fp20;
fp20 = fp52;
_libint2_static_api_inc1_short_(&(stack[((hsi*1+20)*1+lsi)*1]),&(fp20),1);
LIBINT2_REALTYPE fp56;
fp56 = inteval->WP_z[vi] * fp173;
LIBINT2_REALTYPE fp57;
fp57 = inteval->PA_z[vi] * fp159;
LIBINT2_REALTYPE fp55;
fp55 = fp57 + fp56;
LIBINT2_REALTYPE fp19;
fp19 = fp55;
_libint2_static_api_inc1_short_(&(stack[((hsi*1+19)*1+lsi)*1]),&(fp19),1);
LIBINT2_REALTYPE fp181;
fp181 = inteval->WQ_x[vi] * fp156;
LIBINT2_REALTYPE fp182;
fp182 = inteval->QC_x[vi] * fp173;
LIBINT2_REALTYPE fp180;
fp180 = fp182 + fp181;
LIBINT2_REALTYPE fp176;
fp176 = fp180 + fp177;
LIBINT2_REALTYPE fp59;
fp59 = inteval->WP_z[vi] * fp176;
LIBINT2_REALTYPE fp189;
fp189 = inteval->WQ_x[vi] * fp173;
LIBINT2_REALTYPE fp190;
fp190 = inteval->QC_x[vi] * fp159;
LIBINT2_REALTYPE fp188;
fp188 = fp190 + fp189;
LIBINT2_REALTYPE fp183;
fp183 = fp188 + fp184;
LIBINT2_REALTYPE fp60;
fp60 = inteval->PA_z[vi] * fp183;
LIBINT2_REALTYPE fp58;
fp58 = fp60 + fp59;
LIBINT2_REALTYPE fp18;
fp18 = fp58;
_libint2_static_api_inc1_short_(&(stack[((hsi*1+18)*1+lsi)*1]),&(fp18),1);
LIBINT2_REALTYPE fp62;
fp62 = inteval->WP_y[vi] * fp95;
LIBINT2_REALTYPE fp63;
fp63 = inteval->PA_y[vi] * fp99;
LIBINT2_REALTYPE fp61;
fp61 = fp63 + fp62;
LIBINT2_REALTYPE fp17;
fp17 = fp61;
_libint2_static_api_inc1_short_(&(stack[((hsi*1+17)*1+lsi)*1]),&(fp17),1);
LIBINT2_REALTYPE fp73;
fp73 = inteval->WP_y[vi] * fp135;
LIBINT2_REALTYPE fp74;
fp74 = inteval->PA_y[vi] * fp109;
LIBINT2_REALTYPE fp72;
fp72 = fp74 + fp73;
LIBINT2_REALTYPE fp14;
fp14 = fp72;
_libint2_static_api_inc1_short_(&(stack[((hsi*1+16)*1+lsi)*1]),&(fp14),1);
LIBINT2_REALTYPE fp80;
fp80 = inteval->WP_y[vi] * fp115;
LIBINT2_REALTYPE fp81;
fp81 = inteval->PA_y[vi] * fp144;
LIBINT2_REALTYPE fp79;
fp79 = fp81 + fp80;
LIBINT2_REALTYPE fp78;
fp78 = fp79 + fp163;
LIBINT2_REALTYPE fp65;
fp65 = inteval->WQ_z[vi] * fp78;
LIBINT2_REALTYPE fp84;
fp84 = inteval->WP_y[vi] * fp144;
LIBINT2_REALTYPE fp85;
fp85 = inteval->PA_y[vi] * fp129;
LIBINT2_REALTYPE fp83;
fp83 = fp85 + fp84;
LIBINT2_REALTYPE fp82;
fp82 = fp83 + fp168;
LIBINT2_REALTYPE fp12;
fp12 = fp82;
LIBINT2_REALTYPE fp66;
fp66 = inteval->QC_z[vi] * fp12;
LIBINT2_REALTYPE fp64;
fp64 = fp66 + fp65;
LIBINT2_REALTYPE fp16;
fp16 = fp64;
_libint2_static_api_inc1_short_(&(stack[((hsi*1+15)*1+lsi)*1]),&(fp16),1);
LIBINT2_REALTYPE fp68;
fp68 = fp193 * fp144;
LIBINT2_REALTYPE fp70;
fp70 = inteval->WP_y[vi] * fp118;
LIBINT2_REALTYPE fp71;
fp71 = inteval->PA_y[vi] * fp122;
LIBINT2_REALTYPE fp69;
fp69 = fp71 + fp70;
LIBINT2_REALTYPE fp67;
fp67 = fp69 + fp68;
LIBINT2_REALTYPE fp15;
fp15 = fp67;
_libint2_static_api_inc1_short_(&(stack[((hsi*1+14)*1+lsi)*1]),&(fp15),1);
_libint2_static_api_inc1_short_(&(stack[((hsi*1+13)*1+lsi)*1]),&(fp12),1);
LIBINT2_REALTYPE fp207;
fp207 = inteval->WQ_z[vi] * fp156;
LIBINT2_REALTYPE fp208;
fp208 = inteval->QC_z[vi] * fp173;
LIBINT2_REALTYPE fp206;
fp206 = fp208 + fp207;
LIBINT2_REALTYPE fp76;
fp76 = inteval->WP_y[vi] * fp206;
LIBINT2_REALTYPE fp201;
fp201 = inteval->WQ_z[vi] * fp173;
LIBINT2_REALTYPE fp202;
fp202 = inteval->QC_z[vi] * fp159;
LIBINT2_REALTYPE fp200;
fp200 = fp202 + fp201;
LIBINT2_REALTYPE fp77;
fp77 = inteval->PA_y[vi] * fp200;
LIBINT2_REALTYPE fp75;
fp75 = fp77 + fp76;
LIBINT2_REALTYPE fp13;
fp13 = fp75;
_libint2_static_api_inc1_short_(&(stack[((hsi*1+12)*1+lsi)*1]),&(fp13),1);
LIBINT2_REALTYPE fp87;
fp87 = inteval->WQ_x[vi] * fp78;
LIBINT2_REALTYPE fp88;
fp88 = inteval->QC_x[vi] * fp12;
LIBINT2_REALTYPE fp86;
fp86 = fp88 + fp87;
LIBINT2_REALTYPE fp11;
fp11 = fp86;
_libint2_static_api_inc1_short_(&(stack[((hsi*1+11)*1+lsi)*1]),&(fp11),1);
LIBINT2_REALTYPE fp90;
fp90 = inteval->WP_y[vi] * fp173;
LIBINT2_REALTYPE fp91;
fp91 = inteval->PA_y[vi] * fp159;
LIBINT2_REALTYPE fp89;
fp89 = fp91 + fp90;
LIBINT2_REALTYPE fp10;
fp10 = fp89;
_libint2_static_api_inc1_short_(&(stack[((hsi*1+10)*1+lsi)*1]),&(fp10),1);
LIBINT2_REALTYPE fp93;
fp93 = inteval->WP_y[vi] * fp176;
LIBINT2_REALTYPE fp94;
fp94 = inteval->PA_y[vi] * fp183;
LIBINT2_REALTYPE fp92;
fp92 = fp94 + fp93;
LIBINT2_REALTYPE fp9;
fp9 = fp92;
_libint2_static_api_inc1_short_(&(stack[((hsi*1+9)*1+lsi)*1]),&(fp9),1);
LIBINT2_REALTYPE fp104;
fp104 = inteval->WP_x[vi] * fp95;
LIBINT2_REALTYPE fp105;
fp105 = inteval->PA_x[vi] * fp99;
LIBINT2_REALTYPE fp103;
fp103 = fp105 + fp104;
LIBINT2_REALTYPE fp8;
fp8 = fp103;
_libint2_static_api_inc1_short_(&(stack[((hsi*1+8)*1+lsi)*1]),&(fp8),1);
LIBINT2_REALTYPE fp133;
fp133 = inteval->WP_x[vi] * fp135;
LIBINT2_REALTYPE fp134;
fp134 = inteval->PA_x[vi] * fp109;
LIBINT2_REALTYPE fp132;
fp132 = fp134 + fp133;
LIBINT2_REALTYPE fp5;
fp5 = fp132;
_libint2_static_api_inc1_short_(&(stack[((hsi*1+7)*1+lsi)*1]),&(fp5),1);
LIBINT2_REALTYPE fp204;
fp204 = inteval->WQ_z[vi] * fp115;
LIBINT2_REALTYPE fp205;
fp205 = inteval->QC_z[vi] * fp144;
LIBINT2_REALTYPE fp203;
fp203 = fp205 + fp204;
LIBINT2_REALTYPE fp113;
fp113 = inteval->WP_x[vi] * fp203;
LIBINT2_REALTYPE fp198;
fp198 = inteval->WQ_z[vi] * fp144;
LIBINT2_REALTYPE fp199;
fp199 = inteval->QC_z[vi] * fp129;
LIBINT2_REALTYPE fp197;
fp197 = fp199 + fp198;
LIBINT2_REALTYPE fp114;
fp114 = inteval->PA_x[vi] * fp197;
LIBINT2_REALTYPE fp112;
fp112 = fp114 + fp113;
LIBINT2_REALTYPE fp7;
fp7 = fp112;
_libint2_static_api_inc1_short_(&(stack[((hsi*1+6)*1+lsi)*1]),&(fp7),1);
LIBINT2_REALTYPE fp127;
fp127 = inteval->WP_x[vi] * fp118;
LIBINT2_REALTYPE fp128;
fp128 = inteval->PA_x[vi] * fp122;
LIBINT2_REALTYPE fp126;
fp126 = fp128 + fp127;
LIBINT2_REALTYPE fp6;
fp6 = fp126;
_libint2_static_api_inc1_short_(&(stack[((hsi*1+5)*1+lsi)*1]),&(fp6),1);
LIBINT2_REALTYPE fp142;
fp142 = inteval->WP_x[vi] * fp144;
LIBINT2_REALTYPE fp143;
fp143 = inteval->PA_x[vi] * fp129;
LIBINT2_REALTYPE fp141;
fp141 = fp143 + fp142;
LIBINT2_REALTYPE fp3;
fp3 = fp141;
_libint2_static_api_inc1_short_(&(stack[((hsi*1+4)*1+lsi)*1]),&(fp3),1);
LIBINT2_REALTYPE fp165;
fp165 = inteval->WP_x[vi] * fp156;
LIBINT2_REALTYPE fp166;
fp166 = inteval->PA_x[vi] * fp173;
LIBINT2_REALTYPE fp164;
fp164 = fp166 + fp165;
LIBINT2_REALTYPE fp162;
fp162 = fp164 + fp163;
LIBINT2_REALTYPE fp139;
fp139 = inteval->WQ_z[vi] * fp162;
LIBINT2_REALTYPE fp171;
fp171 = inteval->WP_x[vi] * fp173;
LIBINT2_REALTYPE fp172;
fp172 = inteval->PA_x[vi] * fp159;
LIBINT2_REALTYPE fp170;
fp170 = fp172 + fp171;
LIBINT2_REALTYPE fp167;
fp167 = fp170 + fp168;
LIBINT2_REALTYPE fp1;
fp1 = fp167;
LIBINT2_REALTYPE fp140;
fp140 = inteval->QC_z[vi] * fp1;
LIBINT2_REALTYPE fp138;
fp138 = fp140 + fp139;
LIBINT2_REALTYPE fp4;
fp4 = fp138;
_libint2_static_api_inc1_short_(&(stack[((hsi*1+3)*1+lsi)*1]),&(fp4),1);
LIBINT2_REALTYPE fp154;
fp154 = inteval->WQ_y[vi] * fp162;
LIBINT2_REALTYPE fp155;
fp155 = inteval->QC_y[vi] * fp1;
LIBINT2_REALTYPE fp153;
fp153 = fp155 + fp154;
LIBINT2_REALTYPE fp2;
fp2 = fp153;
_libint2_static_api_inc1_short_(&(stack[((hsi*1+2)*1+lsi)*1]),&(fp2),1);
_libint2_static_api_inc1_short_(&(stack[((hsi*1+1)*1+lsi)*1]),&(fp1),1);
LIBINT2_REALTYPE fp192;
fp192 = fp193 * fp173;
LIBINT2_REALTYPE fp195;
fp195 = inteval->WP_x[vi] * fp176;
LIBINT2_REALTYPE fp196;
fp196 = inteval->PA_x[vi] * fp183;
LIBINT2_REALTYPE fp194;
fp194 = fp196 + fp195;
LIBINT2_REALTYPE fp191;
fp191 = fp194 + fp192;
LIBINT2_REALTYPE fp0;
fp0 = fp191;
_libint2_static_api_inc1_short_(&(stack[((hsi*1+0)*1+lsi)*1]),&(fp0),1);
}
}
}
const int hsi = 0;
const int lsi = 0;
const int vi = 0;
/** Number of flops = 209 */
}

#ifdef __cplusplus
};
#endif
