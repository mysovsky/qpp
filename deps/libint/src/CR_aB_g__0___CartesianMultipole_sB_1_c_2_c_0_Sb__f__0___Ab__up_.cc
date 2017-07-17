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

#ifdef __cplusplus
extern "C" {
#endif
void CR_aB_g__0___CartesianMultipole_sB_1_c_2_c_0_Sb__f__0___Ab__up_(const Libint_t* inteval, LIBINT2_REALTYPE* target, const LIBINT2_REALTYPE* src0, const LIBINT2_REALTYPE* src1, const LIBINT2_REALTYPE* src2) {

LIBINT2_REALTYPE*const  stack = target;
{
const int hsi = 0;
{
const int lsi = 0;
{
const int vi = 0;
LIBINT2_REALTYPE fp9;
fp9 = 0 + src2[((hsi*35+16)*1+lsi)*1];
LIBINT2_REALTYPE fp133;
fp133 = 0 + src1[((hsi*35+2)*1+lsi)*1];
LIBINT2_REALTYPE fp132;
fp132 = 0 + src1[((hsi*35+1)*1+lsi)*1];
LIBINT2_REALTYPE fp137;
fp137 = inteval->BO_y[vi] * fp132;
LIBINT2_REALTYPE fp136;
fp136 = fp133 + fp137;
LIBINT2_REALTYPE fp131;
fp131 = 0 + src1[((hsi*35+0)*1+lsi)*1];
LIBINT2_REALTYPE fp135;
fp135 = inteval->BO_y[vi] * fp131;
LIBINT2_REALTYPE fp134;
fp134 = fp132 + fp135;
LIBINT2_REALTYPE fp139;
fp139 = inteval->BO_y[vi] * fp134;
LIBINT2_REALTYPE fp138;
fp138 = fp136 + fp139;
LIBINT2_REALTYPE fp172;
fp172 = 0 + src0[((hsi*35+16)*1+lsi)*1];
LIBINT2_REALTYPE fp169;
fp169 = 0 + src0[((hsi*35+15)*1+lsi)*1];
LIBINT2_REALTYPE fp171;
fp171 = inteval->BO_x[vi] * fp169;
LIBINT2_REALTYPE fp170;
fp170 = fp172 + fp171;
LIBINT2_REALTYPE fp394;
fp394 = fp170 * fp138;
LIBINT2_REALTYPE fp393;
fp393 = fp394 * fp9;
target[((hsi*150+55)*1+lsi)*1] = fp393;
LIBINT2_REALTYPE fp10;
fp10 = 0 + src2[((hsi*35+15)*1+lsi)*1];
LIBINT2_REALTYPE fp126;
fp126 = 0 + src1[((hsi*35+3)*1+lsi)*1];
LIBINT2_REALTYPE fp128;
fp128 = inteval->BO_y[vi] * fp133;
LIBINT2_REALTYPE fp127;
fp127 = fp126 + fp128;
LIBINT2_REALTYPE fp130;
fp130 = inteval->BO_y[vi] * fp136;
LIBINT2_REALTYPE fp129;
fp129 = fp127 + fp130;
LIBINT2_REALTYPE fp396;
fp396 = fp170 * fp129;
LIBINT2_REALTYPE fp395;
fp395 = fp396 * fp10;
target[((hsi*150+54)*1+lsi)*1] = fp395;
LIBINT2_REALTYPE fp11;
fp11 = 0 + src2[((hsi*35+14)*1+lsi)*1];
LIBINT2_REALTYPE fp121;
fp121 = 0 + src1[((hsi*35+4)*1+lsi)*1];
LIBINT2_REALTYPE fp123;
fp123 = inteval->BO_y[vi] * fp126;
LIBINT2_REALTYPE fp122;
fp122 = fp121 + fp123;
LIBINT2_REALTYPE fp125;
fp125 = inteval->BO_y[vi] * fp127;
LIBINT2_REALTYPE fp124;
fp124 = fp122 + fp125;
LIBINT2_REALTYPE fp398;
fp398 = fp170 * fp124;
LIBINT2_REALTYPE fp397;
fp397 = fp398 * fp11;
target[((hsi*150+53)*1+lsi)*1] = fp397;
LIBINT2_REALTYPE fp175;
fp175 = 0 + src0[((hsi*35+17)*1+lsi)*1];
LIBINT2_REALTYPE fp174;
fp174 = inteval->BO_x[vi] * fp172;
LIBINT2_REALTYPE fp173;
fp173 = fp175 + fp174;
LIBINT2_REALTYPE fp400;
fp400 = fp173 * fp138;
LIBINT2_REALTYPE fp399;
fp399 = fp400 * fp10;
target[((hsi*150+52)*1+lsi)*1] = fp399;
LIBINT2_REALTYPE fp402;
fp402 = fp173 * fp129;
LIBINT2_REALTYPE fp401;
fp401 = fp402 * fp11;
target[((hsi*150+51)*1+lsi)*1] = fp401;
LIBINT2_REALTYPE fp178;
fp178 = inteval->BO_x[vi] * fp175;
LIBINT2_REALTYPE fp176;
fp176 = 0 + src0[((hsi*35+18)*1+lsi)*1];
LIBINT2_REALTYPE fp177;
fp177 = fp176 + fp178;
LIBINT2_REALTYPE fp404;
fp404 = fp177 * fp138;
LIBINT2_REALTYPE fp403;
fp403 = fp404 * fp11;
target[((hsi*150+50)*1+lsi)*1] = fp403;
LIBINT2_REALTYPE fp12;
fp12 = 0 + src2[((hsi*35+10)*1+lsi)*1];
LIBINT2_REALTYPE fp109;
fp109 = 0 + src1[((hsi*35+9)*1+lsi)*1];
LIBINT2_REALTYPE fp108;
fp108 = 0 + src1[((hsi*35+8)*1+lsi)*1];
LIBINT2_REALTYPE fp113;
fp113 = inteval->BO_y[vi] * fp108;
LIBINT2_REALTYPE fp112;
fp112 = fp109 + fp113;
LIBINT2_REALTYPE fp107;
fp107 = 0 + src1[((hsi*35+7)*1+lsi)*1];
LIBINT2_REALTYPE fp111;
fp111 = inteval->BO_y[vi] * fp107;
LIBINT2_REALTYPE fp110;
fp110 = fp108 + fp111;
LIBINT2_REALTYPE fp115;
fp115 = inteval->BO_y[vi] * fp110;
LIBINT2_REALTYPE fp114;
fp114 = fp112 + fp115;
LIBINT2_REALTYPE fp166;
fp166 = 0 + src0[((hsi*35+14)*1+lsi)*1];
LIBINT2_REALTYPE fp168;
fp168 = inteval->BO_x[vi] * fp166;
LIBINT2_REALTYPE fp167;
fp167 = fp169 + fp168;
LIBINT2_REALTYPE fp406;
fp406 = fp167 * fp114;
LIBINT2_REALTYPE fp405;
fp405 = fp406 * fp12;
target[((hsi*150+49)*1+lsi)*1] = fp405;
LIBINT2_REALTYPE fp13;
fp13 = 0 + src2[((hsi*35+9)*1+lsi)*1];
LIBINT2_REALTYPE fp102;
fp102 = 0 + src1[((hsi*35+10)*1+lsi)*1];
LIBINT2_REALTYPE fp104;
fp104 = inteval->BO_y[vi] * fp109;
LIBINT2_REALTYPE fp103;
fp103 = fp102 + fp104;
LIBINT2_REALTYPE fp106;
fp106 = inteval->BO_y[vi] * fp112;
LIBINT2_REALTYPE fp105;
fp105 = fp103 + fp106;
LIBINT2_REALTYPE fp408;
fp408 = fp167 * fp105;
LIBINT2_REALTYPE fp407;
fp407 = fp408 * fp13;
target[((hsi*150+48)*1+lsi)*1] = fp407;
LIBINT2_REALTYPE fp19;
fp19 = 0 + src2[((hsi*35+0)*1+lsi)*1];
LIBINT2_REALTYPE fp201;
fp201 = 0 + src0[((hsi*35+31)*1+lsi)*1];
LIBINT2_REALTYPE fp204;
fp204 = inteval->BO_x[vi] * fp201;
LIBINT2_REALTYPE fp202;
fp202 = 0 + src0[((hsi*35+32)*1+lsi)*1];
LIBINT2_REALTYPE fp203;
fp203 = fp202 + fp204;
LIBINT2_REALTYPE fp504;
fp504 = fp203 * fp138;
LIBINT2_REALTYPE fp503;
fp503 = fp504 * fp19;
target[((hsi*150+0)*1+lsi)*1] = fp503;
LIBINT2_REALTYPE fp15;
fp15 = 0 + src2[((hsi*35+7)*1+lsi)*1];
LIBINT2_REALTYPE fp97;
fp97 = 0 + src1[((hsi*35+11)*1+lsi)*1];
LIBINT2_REALTYPE fp99;
fp99 = inteval->BO_y[vi] * fp102;
LIBINT2_REALTYPE fp98;
fp98 = fp97 + fp99;
LIBINT2_REALTYPE fp96;
fp96 = inteval->BO_y[vi] * fp98;
LIBINT2_REALTYPE fp94;
fp94 = inteval->BO_y[vi] * fp97;
LIBINT2_REALTYPE fp92;
fp92 = 0 + src1[((hsi*35+12)*1+lsi)*1];
LIBINT2_REALTYPE fp93;
fp93 = fp92 + fp94;
LIBINT2_REALTYPE fp95;
fp95 = fp93 + fp96;
LIBINT2_REALTYPE fp412;
fp412 = fp167 * fp95;
LIBINT2_REALTYPE fp411;
fp411 = fp412 * fp15;
target[((hsi*150+46)*1+lsi)*1] = fp411;
LIBINT2_REALTYPE fp414;
fp414 = fp170 * fp114;
LIBINT2_REALTYPE fp413;
fp413 = fp414 * fp13;
target[((hsi*150+45)*1+lsi)*1] = fp413;
LIBINT2_REALTYPE fp14;
fp14 = 0 + src2[((hsi*35+8)*1+lsi)*1];
LIBINT2_REALTYPE fp416;
fp416 = fp170 * fp105;
LIBINT2_REALTYPE fp415;
fp415 = fp416 * fp14;
target[((hsi*150+44)*1+lsi)*1] = fp415;
LIBINT2_REALTYPE fp101;
fp101 = inteval->BO_y[vi] * fp103;
LIBINT2_REALTYPE fp100;
fp100 = fp98 + fp101;
LIBINT2_REALTYPE fp418;
fp418 = fp170 * fp100;
LIBINT2_REALTYPE fp417;
fp417 = fp418 * fp15;
target[((hsi*150+43)*1+lsi)*1] = fp417;
LIBINT2_REALTYPE fp420;
fp420 = fp173 * fp114;
LIBINT2_REALTYPE fp419;
fp419 = fp420 * fp14;
target[((hsi*150+42)*1+lsi)*1] = fp419;
LIBINT2_REALTYPE fp422;
fp422 = fp173 * fp105;
LIBINT2_REALTYPE fp421;
fp421 = fp422 * fp15;
target[((hsi*150+41)*1+lsi)*1] = fp421;
LIBINT2_REALTYPE fp424;
fp424 = fp177 * fp114;
LIBINT2_REALTYPE fp423;
fp423 = fp424 * fp15;
target[((hsi*150+40)*1+lsi)*1] = fp423;
LIBINT2_REALTYPE fp16;
fp16 = 0 + src2[((hsi*35+3)*1+lsi)*1];
LIBINT2_REALTYPE fp85;
fp85 = 0 + src1[((hsi*35+16)*1+lsi)*1];
LIBINT2_REALTYPE fp84;
fp84 = 0 + src1[((hsi*35+15)*1+lsi)*1];
LIBINT2_REALTYPE fp89;
fp89 = inteval->BO_y[vi] * fp84;
LIBINT2_REALTYPE fp88;
fp88 = fp85 + fp89;
LIBINT2_REALTYPE fp83;
fp83 = 0 + src1[((hsi*35+14)*1+lsi)*1];
LIBINT2_REALTYPE fp87;
fp87 = inteval->BO_y[vi] * fp83;
LIBINT2_REALTYPE fp86;
fp86 = fp84 + fp87;
LIBINT2_REALTYPE fp91;
fp91 = inteval->BO_y[vi] * fp86;
LIBINT2_REALTYPE fp90;
fp90 = fp88 + fp91;
LIBINT2_REALTYPE fp426;
fp426 = fp167 * fp90;
LIBINT2_REALTYPE fp425;
fp425 = fp426 * fp16;
target[((hsi*150+39)*1+lsi)*1] = fp425;
LIBINT2_REALTYPE fp410;
fp410 = fp167 * fp100;
LIBINT2_REALTYPE fp409;
fp409 = fp410 * fp14;
target[((hsi*150+47)*1+lsi)*1] = fp409;
LIBINT2_REALTYPE fp159;
fp159 = 0 + src0[((hsi*35+9)*1+lsi)*1];
LIBINT2_REALTYPE fp156;
fp156 = 0 + src0[((hsi*35+8)*1+lsi)*1];
LIBINT2_REALTYPE fp158;
fp158 = inteval->BO_x[vi] * fp156;
LIBINT2_REALTYPE fp157;
fp157 = fp159 + fp158;
LIBINT2_REALTYPE fp78;
fp78 = 0 + src1[((hsi*35+17)*1+lsi)*1];
LIBINT2_REALTYPE fp80;
fp80 = inteval->BO_y[vi] * fp85;
LIBINT2_REALTYPE fp79;
fp79 = fp78 + fp80;
LIBINT2_REALTYPE fp82;
fp82 = inteval->BO_y[vi] * fp88;
LIBINT2_REALTYPE fp81;
fp81 = fp79 + fp82;
LIBINT2_REALTYPE fp356;
fp356 = fp157 * fp81;
LIBINT2_REALTYPE fp355;
fp355 = fp356 * fp14;
target[((hsi*150+74)*1+lsi)*1] = fp355;
LIBINT2_REALTYPE fp73;
fp73 = 0 + src1[((hsi*35+18)*1+lsi)*1];
LIBINT2_REALTYPE fp75;
fp75 = inteval->BO_y[vi] * fp78;
LIBINT2_REALTYPE fp74;
fp74 = fp73 + fp75;
LIBINT2_REALTYPE fp77;
fp77 = inteval->BO_y[vi] * fp79;
LIBINT2_REALTYPE fp76;
fp76 = fp74 + fp77;
LIBINT2_REALTYPE fp358;
fp358 = fp157 * fp76;
LIBINT2_REALTYPE fp357;
fp357 = fp358 * fp15;
target[((hsi*150+73)*1+lsi)*1] = fp357;
LIBINT2_REALTYPE fp162;
fp162 = 0 + src0[((hsi*35+10)*1+lsi)*1];
LIBINT2_REALTYPE fp161;
fp161 = inteval->BO_x[vi] * fp159;
LIBINT2_REALTYPE fp160;
fp160 = fp162 + fp161;
LIBINT2_REALTYPE fp360;
fp360 = fp160 * fp90;
LIBINT2_REALTYPE fp359;
fp359 = fp360 * fp14;
target[((hsi*150+72)*1+lsi)*1] = fp359;
LIBINT2_REALTYPE fp362;
fp362 = fp160 * fp81;
LIBINT2_REALTYPE fp361;
fp361 = fp362 * fp15;
target[((hsi*150+71)*1+lsi)*1] = fp361;
LIBINT2_REALTYPE fp165;
fp165 = inteval->BO_x[vi] * fp162;
LIBINT2_REALTYPE fp163;
fp163 = 0 + src0[((hsi*35+11)*1+lsi)*1];
LIBINT2_REALTYPE fp164;
fp164 = fp163 + fp165;
LIBINT2_REALTYPE fp364;
fp364 = fp164 * fp90;
LIBINT2_REALTYPE fp363;
fp363 = fp364 * fp15;
target[((hsi*150+70)*1+lsi)*1] = fp363;
LIBINT2_REALTYPE fp153;
fp153 = 0 + src0[((hsi*35+7)*1+lsi)*1];
LIBINT2_REALTYPE fp155;
fp155 = inteval->BO_x[vi] * fp153;
LIBINT2_REALTYPE fp154;
fp154 = fp156 + fp155;
LIBINT2_REALTYPE fp61;
fp61 = 0 + src1[((hsi*35+23)*1+lsi)*1];
LIBINT2_REALTYPE fp60;
fp60 = 0 + src1[((hsi*35+22)*1+lsi)*1];
LIBINT2_REALTYPE fp65;
fp65 = inteval->BO_y[vi] * fp60;
LIBINT2_REALTYPE fp64;
fp64 = fp61 + fp65;
LIBINT2_REALTYPE fp59;
fp59 = 0 + src1[((hsi*35+21)*1+lsi)*1];
LIBINT2_REALTYPE fp63;
fp63 = inteval->BO_y[vi] * fp59;
LIBINT2_REALTYPE fp62;
fp62 = fp60 + fp63;
LIBINT2_REALTYPE fp67;
fp67 = inteval->BO_y[vi] * fp62;
LIBINT2_REALTYPE fp66;
fp66 = fp64 + fp67;
LIBINT2_REALTYPE fp366;
fp366 = fp154 * fp66;
LIBINT2_REALTYPE fp365;
fp365 = fp366 * fp16;
target[((hsi*150+69)*1+lsi)*1] = fp365;
LIBINT2_REALTYPE fp17;
fp17 = 0 + src2[((hsi*35+2)*1+lsi)*1];
LIBINT2_REALTYPE fp54;
fp54 = 0 + src1[((hsi*35+24)*1+lsi)*1];
LIBINT2_REALTYPE fp56;
fp56 = inteval->BO_y[vi] * fp61;
LIBINT2_REALTYPE fp55;
fp55 = fp54 + fp56;
LIBINT2_REALTYPE fp58;
fp58 = inteval->BO_y[vi] * fp64;
LIBINT2_REALTYPE fp57;
fp57 = fp55 + fp58;
LIBINT2_REALTYPE fp368;
fp368 = fp154 * fp57;
LIBINT2_REALTYPE fp367;
fp367 = fp368 * fp17;
target[((hsi*150+68)*1+lsi)*1] = fp367;
LIBINT2_REALTYPE fp18;
fp18 = 0 + src2[((hsi*35+1)*1+lsi)*1];
LIBINT2_REALTYPE fp49;
fp49 = 0 + src1[((hsi*35+25)*1+lsi)*1];
LIBINT2_REALTYPE fp51;
fp51 = inteval->BO_y[vi] * fp54;
LIBINT2_REALTYPE fp50;
fp50 = fp49 + fp51;
LIBINT2_REALTYPE fp53;
fp53 = inteval->BO_y[vi] * fp55;
LIBINT2_REALTYPE fp52;
fp52 = fp50 + fp53;
LIBINT2_REALTYPE fp370;
fp370 = fp154 * fp52;
LIBINT2_REALTYPE fp369;
fp369 = fp370 * fp18;
target[((hsi*150+67)*1+lsi)*1] = fp369;
LIBINT2_REALTYPE fp120;
fp120 = inteval->BO_y[vi] * fp122;
LIBINT2_REALTYPE fp118;
fp118 = inteval->BO_y[vi] * fp121;
LIBINT2_REALTYPE fp116;
fp116 = 0 + src1[((hsi*35+5)*1+lsi)*1];
LIBINT2_REALTYPE fp117;
fp117 = fp116 + fp118;
LIBINT2_REALTYPE fp119;
fp119 = fp117 + fp120;
LIBINT2_REALTYPE fp392;
fp392 = fp167 * fp119;
LIBINT2_REALTYPE fp391;
fp391 = fp392 * fp11;
target[((hsi*150+56)*1+lsi)*1] = fp391;
LIBINT2_REALTYPE fp374;
fp374 = fp157 * fp66;
LIBINT2_REALTYPE fp373;
fp373 = fp374 * fp17;
target[((hsi*150+65)*1+lsi)*1] = fp373;
LIBINT2_REALTYPE fp376;
fp376 = fp157 * fp57;
LIBINT2_REALTYPE fp375;
fp375 = fp376 * fp18;
target[((hsi*150+64)*1+lsi)*1] = fp375;
LIBINT2_REALTYPE fp378;
fp378 = fp157 * fp52;
LIBINT2_REALTYPE fp377;
fp377 = fp378 * fp19;
target[((hsi*150+63)*1+lsi)*1] = fp377;
LIBINT2_REALTYPE fp380;
fp380 = fp160 * fp66;
LIBINT2_REALTYPE fp379;
fp379 = fp380 * fp18;
target[((hsi*150+62)*1+lsi)*1] = fp379;
LIBINT2_REALTYPE fp382;
fp382 = fp160 * fp57;
LIBINT2_REALTYPE fp381;
fp381 = fp382 * fp19;
target[((hsi*150+61)*1+lsi)*1] = fp381;
LIBINT2_REALTYPE fp384;
fp384 = fp164 * fp66;
LIBINT2_REALTYPE fp383;
fp383 = fp384 * fp19;
target[((hsi*150+60)*1+lsi)*1] = fp383;
LIBINT2_REALTYPE fp8;
fp8 = 0 + src2[((hsi*35+17)*1+lsi)*1];
LIBINT2_REALTYPE fp386;
fp386 = fp167 * fp138;
LIBINT2_REALTYPE fp385;
fp385 = fp386 * fp8;
target[((hsi*150+59)*1+lsi)*1] = fp385;
LIBINT2_REALTYPE fp388;
fp388 = fp167 * fp129;
LIBINT2_REALTYPE fp387;
fp387 = fp388 * fp9;
target[((hsi*150+58)*1+lsi)*1] = fp387;
LIBINT2_REALTYPE fp390;
fp390 = fp167 * fp124;
LIBINT2_REALTYPE fp389;
fp389 = fp390 * fp10;
target[((hsi*150+57)*1+lsi)*1] = fp389;
LIBINT2_REALTYPE fp48;
fp48 = inteval->BO_y[vi] * fp50;
LIBINT2_REALTYPE fp46;
fp46 = inteval->BO_y[vi] * fp49;
LIBINT2_REALTYPE fp44;
fp44 = 0 + src1[((hsi*35+26)*1+lsi)*1];
LIBINT2_REALTYPE fp45;
fp45 = fp44 + fp46;
LIBINT2_REALTYPE fp47;
fp47 = fp45 + fp48;
LIBINT2_REALTYPE fp372;
fp372 = fp154 * fp47;
LIBINT2_REALTYPE fp371;
fp371 = fp372 * fp19;
target[((hsi*150+66)*1+lsi)*1] = fp371;
LIBINT2_REALTYPE fp182;
fp182 = 0 + src0[((hsi*35+22)*1+lsi)*1];
LIBINT2_REALTYPE fp179;
fp179 = 0 + src0[((hsi*35+21)*1+lsi)*1];
LIBINT2_REALTYPE fp181;
fp181 = inteval->BO_x[vi] * fp179;
LIBINT2_REALTYPE fp180;
fp180 = fp182 + fp181;
LIBINT2_REALTYPE fp468;
fp468 = fp180 * fp105;
LIBINT2_REALTYPE fp467;
fp467 = fp468 * fp17;
target[((hsi*150+18)*1+lsi)*1] = fp467;
LIBINT2_REALTYPE fp470;
fp470 = fp180 * fp100;
LIBINT2_REALTYPE fp469;
fp469 = fp470 * fp18;
target[((hsi*150+17)*1+lsi)*1] = fp469;
LIBINT2_REALTYPE fp472;
fp472 = fp180 * fp95;
LIBINT2_REALTYPE fp471;
fp471 = fp472 * fp19;
target[((hsi*150+16)*1+lsi)*1] = fp471;
LIBINT2_REALTYPE fp185;
fp185 = 0 + src0[((hsi*35+23)*1+lsi)*1];
LIBINT2_REALTYPE fp184;
fp184 = inteval->BO_x[vi] * fp182;
LIBINT2_REALTYPE fp183;
fp183 = fp185 + fp184;
LIBINT2_REALTYPE fp474;
fp474 = fp183 * fp114;
LIBINT2_REALTYPE fp473;
fp473 = fp474 * fp17;
target[((hsi*150+15)*1+lsi)*1] = fp473;
LIBINT2_REALTYPE fp476;
fp476 = fp183 * fp105;
LIBINT2_REALTYPE fp475;
fp475 = fp476 * fp18;
target[((hsi*150+14)*1+lsi)*1] = fp475;
LIBINT2_REALTYPE fp478;
fp478 = fp183 * fp100;
LIBINT2_REALTYPE fp477;
fp477 = fp478 * fp19;
target[((hsi*150+13)*1+lsi)*1] = fp477;
LIBINT2_REALTYPE fp188;
fp188 = 0 + src0[((hsi*35+24)*1+lsi)*1];
LIBINT2_REALTYPE fp187;
fp187 = inteval->BO_x[vi] * fp185;
LIBINT2_REALTYPE fp186;
fp186 = fp188 + fp187;
LIBINT2_REALTYPE fp480;
fp480 = fp186 * fp114;
LIBINT2_REALTYPE fp479;
fp479 = fp480 * fp18;
target[((hsi*150+12)*1+lsi)*1] = fp479;
LIBINT2_REALTYPE fp482;
fp482 = fp186 * fp105;
LIBINT2_REALTYPE fp481;
fp481 = fp482 * fp19;
target[((hsi*150+11)*1+lsi)*1] = fp481;
LIBINT2_REALTYPE fp428;
fp428 = fp167 * fp81;
LIBINT2_REALTYPE fp427;
fp427 = fp428 * fp17;
target[((hsi*150+38)*1+lsi)*1] = fp427;
LIBINT2_REALTYPE fp195;
fp195 = 0 + src0[((hsi*35+29)*1+lsi)*1];
LIBINT2_REALTYPE fp192;
fp192 = 0 + src0[((hsi*35+28)*1+lsi)*1];
LIBINT2_REALTYPE fp194;
fp194 = inteval->BO_x[vi] * fp192;
LIBINT2_REALTYPE fp193;
fp193 = fp195 + fp194;
LIBINT2_REALTYPE fp486;
fp486 = fp193 * fp138;
LIBINT2_REALTYPE fp485;
fp485 = fp486 * fp16;
target[((hsi*150+9)*1+lsi)*1] = fp485;
LIBINT2_REALTYPE fp488;
fp488 = fp193 * fp129;
LIBINT2_REALTYPE fp487;
fp487 = fp488 * fp17;
target[((hsi*150+8)*1+lsi)*1] = fp487;
LIBINT2_REALTYPE fp490;
fp490 = fp193 * fp124;
LIBINT2_REALTYPE fp489;
fp489 = fp490 * fp18;
target[((hsi*150+7)*1+lsi)*1] = fp489;
LIBINT2_REALTYPE fp492;
fp492 = fp193 * fp119;
LIBINT2_REALTYPE fp491;
fp491 = fp492 * fp19;
target[((hsi*150+6)*1+lsi)*1] = fp491;
LIBINT2_REALTYPE fp198;
fp198 = 0 + src0[((hsi*35+30)*1+lsi)*1];
LIBINT2_REALTYPE fp197;
fp197 = inteval->BO_x[vi] * fp195;
LIBINT2_REALTYPE fp196;
fp196 = fp198 + fp197;
LIBINT2_REALTYPE fp494;
fp494 = fp196 * fp138;
LIBINT2_REALTYPE fp493;
fp493 = fp494 * fp17;
target[((hsi*150+5)*1+lsi)*1] = fp493;
LIBINT2_REALTYPE fp496;
fp496 = fp196 * fp129;
LIBINT2_REALTYPE fp495;
fp495 = fp496 * fp18;
target[((hsi*150+4)*1+lsi)*1] = fp495;
LIBINT2_REALTYPE fp498;
fp498 = fp196 * fp124;
LIBINT2_REALTYPE fp497;
fp497 = fp498 * fp19;
target[((hsi*150+3)*1+lsi)*1] = fp497;
LIBINT2_REALTYPE fp200;
fp200 = inteval->BO_x[vi] * fp198;
LIBINT2_REALTYPE fp199;
fp199 = fp201 + fp200;
LIBINT2_REALTYPE fp500;
fp500 = fp199 * fp138;
LIBINT2_REALTYPE fp499;
fp499 = fp500 * fp18;
target[((hsi*150+2)*1+lsi)*1] = fp499;
LIBINT2_REALTYPE fp502;
fp502 = fp199 * fp129;
LIBINT2_REALTYPE fp501;
fp501 = fp502 * fp19;
target[((hsi*150+1)*1+lsi)*1] = fp501;
LIBINT2_REALTYPE fp191;
fp191 = inteval->BO_x[vi] * fp188;
LIBINT2_REALTYPE fp189;
fp189 = 0 + src0[((hsi*35+25)*1+lsi)*1];
LIBINT2_REALTYPE fp190;
fp190 = fp189 + fp191;
LIBINT2_REALTYPE fp484;
fp484 = fp190 * fp114;
LIBINT2_REALTYPE fp483;
fp483 = fp484 * fp19;
target[((hsi*150+10)*1+lsi)*1] = fp483;
LIBINT2_REALTYPE fp430;
fp430 = fp167 * fp76;
LIBINT2_REALTYPE fp429;
fp429 = fp430 * fp18;
target[((hsi*150+37)*1+lsi)*1] = fp429;
LIBINT2_REALTYPE fp72;
fp72 = inteval->BO_y[vi] * fp74;
LIBINT2_REALTYPE fp70;
fp70 = inteval->BO_y[vi] * fp73;
LIBINT2_REALTYPE fp68;
fp68 = 0 + src1[((hsi*35+19)*1+lsi)*1];
LIBINT2_REALTYPE fp69;
fp69 = fp68 + fp70;
LIBINT2_REALTYPE fp71;
fp71 = fp69 + fp72;
LIBINT2_REALTYPE fp432;
fp432 = fp167 * fp71;
LIBINT2_REALTYPE fp431;
fp431 = fp432 * fp19;
target[((hsi*150+36)*1+lsi)*1] = fp431;
LIBINT2_REALTYPE fp434;
fp434 = fp170 * fp90;
LIBINT2_REALTYPE fp433;
fp433 = fp434 * fp17;
target[((hsi*150+35)*1+lsi)*1] = fp433;
LIBINT2_REALTYPE fp436;
fp436 = fp170 * fp81;
LIBINT2_REALTYPE fp435;
fp435 = fp436 * fp18;
target[((hsi*150+34)*1+lsi)*1] = fp435;
LIBINT2_REALTYPE fp438;
fp438 = fp170 * fp76;
LIBINT2_REALTYPE fp437;
fp437 = fp438 * fp19;
target[((hsi*150+33)*1+lsi)*1] = fp437;
LIBINT2_REALTYPE fp440;
fp440 = fp173 * fp90;
LIBINT2_REALTYPE fp439;
fp439 = fp440 * fp18;
target[((hsi*150+32)*1+lsi)*1] = fp439;
LIBINT2_REALTYPE fp442;
fp442 = fp173 * fp81;
LIBINT2_REALTYPE fp441;
fp441 = fp442 * fp19;
target[((hsi*150+31)*1+lsi)*1] = fp441;
LIBINT2_REALTYPE fp444;
fp444 = fp177 * fp90;
LIBINT2_REALTYPE fp443;
fp443 = fp444 * fp19;
target[((hsi*150+30)*1+lsi)*1] = fp443;
LIBINT2_REALTYPE fp466;
fp466 = fp180 * fp114;
LIBINT2_REALTYPE fp465;
fp465 = fp466 * fp16;
target[((hsi*150+19)*1+lsi)*1] = fp465;
LIBINT2_REALTYPE fp448;
fp448 = fp180 * fp129;
LIBINT2_REALTYPE fp447;
fp447 = fp448 * fp13;
target[((hsi*150+28)*1+lsi)*1] = fp447;
LIBINT2_REALTYPE fp450;
fp450 = fp180 * fp124;
LIBINT2_REALTYPE fp449;
fp449 = fp450 * fp14;
target[((hsi*150+27)*1+lsi)*1] = fp449;
LIBINT2_REALTYPE fp452;
fp452 = fp180 * fp119;
LIBINT2_REALTYPE fp451;
fp451 = fp452 * fp15;
target[((hsi*150+26)*1+lsi)*1] = fp451;
LIBINT2_REALTYPE fp454;
fp454 = fp183 * fp138;
LIBINT2_REALTYPE fp453;
fp453 = fp454 * fp13;
target[((hsi*150+25)*1+lsi)*1] = fp453;
LIBINT2_REALTYPE fp456;
fp456 = fp183 * fp129;
LIBINT2_REALTYPE fp455;
fp455 = fp456 * fp14;
target[((hsi*150+24)*1+lsi)*1] = fp455;
LIBINT2_REALTYPE fp458;
fp458 = fp183 * fp124;
LIBINT2_REALTYPE fp457;
fp457 = fp458 * fp15;
target[((hsi*150+23)*1+lsi)*1] = fp457;
LIBINT2_REALTYPE fp460;
fp460 = fp186 * fp138;
LIBINT2_REALTYPE fp459;
fp459 = fp460 * fp14;
target[((hsi*150+22)*1+lsi)*1] = fp459;
LIBINT2_REALTYPE fp462;
fp462 = fp186 * fp129;
LIBINT2_REALTYPE fp461;
fp461 = fp462 * fp15;
target[((hsi*150+21)*1+lsi)*1] = fp461;
LIBINT2_REALTYPE fp464;
fp464 = fp190 * fp138;
LIBINT2_REALTYPE fp463;
fp463 = fp464 * fp15;
target[((hsi*150+20)*1+lsi)*1] = fp463;
LIBINT2_REALTYPE fp446;
fp446 = fp180 * fp138;
LIBINT2_REALTYPE fp445;
fp445 = fp446 * fp12;
target[((hsi*150+29)*1+lsi)*1] = fp445;
LIBINT2_REALTYPE fp7;
fp7 = 0 + src2[((hsi*35+21)*1+lsi)*1];
LIBINT2_REALTYPE fp149;
fp149 = 0 + src0[((hsi*35+3)*1+lsi)*1];
LIBINT2_REALTYPE fp152;
fp152 = inteval->BO_x[vi] * fp149;
LIBINT2_REALTYPE fp150;
fp150 = 0 + src0[((hsi*35+4)*1+lsi)*1];
LIBINT2_REALTYPE fp151;
fp151 = fp150 + fp152;
LIBINT2_REALTYPE fp244;
fp244 = fp151 * fp114;
LIBINT2_REALTYPE fp243;
fp243 = fp244 * fp7;
target[((hsi*150+130)*1+lsi)*1] = fp243;
LIBINT2_REALTYPE fp143;
fp143 = 0 + src0[((hsi*35+1)*1+lsi)*1];
LIBINT2_REALTYPE fp140;
fp140 = 0 + src0[((hsi*35+0)*1+lsi)*1];
LIBINT2_REALTYPE fp142;
fp142 = inteval->BO_x[vi] * fp140;
LIBINT2_REALTYPE fp141;
fp141 = fp143 + fp142;
LIBINT2_REALTYPE fp246;
fp246 = fp141 * fp90;
LIBINT2_REALTYPE fp245;
fp245 = fp246 * fp8;
target[((hsi*150+129)*1+lsi)*1] = fp245;
LIBINT2_REALTYPE fp248;
fp248 = fp141 * fp81;
LIBINT2_REALTYPE fp247;
fp247 = fp248 * fp9;
target[((hsi*150+128)*1+lsi)*1] = fp247;
LIBINT2_REALTYPE fp250;
fp250 = fp141 * fp76;
LIBINT2_REALTYPE fp249;
fp249 = fp250 * fp10;
target[((hsi*150+127)*1+lsi)*1] = fp249;
LIBINT2_REALTYPE fp252;
fp252 = fp141 * fp71;
LIBINT2_REALTYPE fp251;
fp251 = fp252 * fp11;
target[((hsi*150+126)*1+lsi)*1] = fp251;
LIBINT2_REALTYPE fp146;
fp146 = 0 + src0[((hsi*35+2)*1+lsi)*1];
LIBINT2_REALTYPE fp145;
fp145 = inteval->BO_x[vi] * fp143;
LIBINT2_REALTYPE fp144;
fp144 = fp146 + fp145;
LIBINT2_REALTYPE fp254;
fp254 = fp144 * fp90;
LIBINT2_REALTYPE fp253;
fp253 = fp254 * fp9;
target[((hsi*150+125)*1+lsi)*1] = fp253;
LIBINT2_REALTYPE fp256;
fp256 = fp144 * fp81;
LIBINT2_REALTYPE fp255;
fp255 = fp256 * fp10;
target[((hsi*150+124)*1+lsi)*1] = fp255;
LIBINT2_REALTYPE fp258;
fp258 = fp144 * fp76;
LIBINT2_REALTYPE fp257;
fp257 = fp258 * fp11;
target[((hsi*150+123)*1+lsi)*1] = fp257;
LIBINT2_REALTYPE fp354;
fp354 = fp157 * fp90;
LIBINT2_REALTYPE fp353;
fp353 = fp354 * fp13;
target[((hsi*150+75)*1+lsi)*1] = fp353;
LIBINT2_REALTYPE fp148;
fp148 = inteval->BO_x[vi] * fp146;
LIBINT2_REALTYPE fp147;
fp147 = fp149 + fp148;
LIBINT2_REALTYPE fp262;
fp262 = fp147 * fp81;
LIBINT2_REALTYPE fp261;
fp261 = fp262 * fp11;
target[((hsi*150+121)*1+lsi)*1] = fp261;
LIBINT2_REALTYPE fp264;
fp264 = fp151 * fp90;
LIBINT2_REALTYPE fp263;
fp263 = fp264 * fp11;
target[((hsi*150+120)*1+lsi)*1] = fp263;
LIBINT2_REALTYPE fp266;
fp266 = fp141 * fp66;
LIBINT2_REALTYPE fp265;
fp265 = fp266 * fp12;
target[((hsi*150+119)*1+lsi)*1] = fp265;
LIBINT2_REALTYPE fp268;
fp268 = fp141 * fp57;
LIBINT2_REALTYPE fp267;
fp267 = fp268 * fp13;
target[((hsi*150+118)*1+lsi)*1] = fp267;
LIBINT2_REALTYPE fp270;
fp270 = fp141 * fp52;
LIBINT2_REALTYPE fp269;
fp269 = fp270 * fp14;
target[((hsi*150+117)*1+lsi)*1] = fp269;
LIBINT2_REALTYPE fp272;
fp272 = fp141 * fp47;
LIBINT2_REALTYPE fp271;
fp271 = fp272 * fp15;
target[((hsi*150+116)*1+lsi)*1] = fp271;
LIBINT2_REALTYPE fp274;
fp274 = fp144 * fp66;
LIBINT2_REALTYPE fp273;
fp273 = fp274 * fp13;
target[((hsi*150+115)*1+lsi)*1] = fp273;
LIBINT2_REALTYPE fp276;
fp276 = fp144 * fp57;
LIBINT2_REALTYPE fp275;
fp275 = fp276 * fp14;
target[((hsi*150+114)*1+lsi)*1] = fp275;
LIBINT2_REALTYPE fp260;
fp260 = fp147 * fp90;
LIBINT2_REALTYPE fp259;
fp259 = fp260 * fp10;
target[((hsi*150+122)*1+lsi)*1] = fp259;
LIBINT2_REALTYPE fp0;
fp0 = 0 + src2[((hsi*35+31)*1+lsi)*1];
LIBINT2_REALTYPE fp206;
fp206 = fp141 * fp138;
LIBINT2_REALTYPE fp205;
fp205 = fp206 * fp0;
target[((hsi*150+149)*1+lsi)*1] = fp205;
LIBINT2_REALTYPE fp1;
fp1 = 0 + src2[((hsi*35+30)*1+lsi)*1];
LIBINT2_REALTYPE fp208;
fp208 = fp141 * fp129;
LIBINT2_REALTYPE fp207;
fp207 = fp208 * fp1;
target[((hsi*150+148)*1+lsi)*1] = fp207;
LIBINT2_REALTYPE fp2;
fp2 = 0 + src2[((hsi*35+29)*1+lsi)*1];
LIBINT2_REALTYPE fp210;
fp210 = fp141 * fp124;
LIBINT2_REALTYPE fp209;
fp209 = fp210 * fp2;
target[((hsi*150+147)*1+lsi)*1] = fp209;
LIBINT2_REALTYPE fp3;
fp3 = 0 + src2[((hsi*35+28)*1+lsi)*1];
LIBINT2_REALTYPE fp212;
fp212 = fp141 * fp119;
LIBINT2_REALTYPE fp211;
fp211 = fp212 * fp3;
target[((hsi*150+146)*1+lsi)*1] = fp211;
LIBINT2_REALTYPE fp214;
fp214 = fp144 * fp138;
LIBINT2_REALTYPE fp213;
fp213 = fp214 * fp1;
target[((hsi*150+145)*1+lsi)*1] = fp213;
LIBINT2_REALTYPE fp216;
fp216 = fp144 * fp129;
LIBINT2_REALTYPE fp215;
fp215 = fp216 * fp2;
target[((hsi*150+144)*1+lsi)*1] = fp215;
LIBINT2_REALTYPE fp218;
fp218 = fp144 * fp124;
LIBINT2_REALTYPE fp217;
fp217 = fp218 * fp3;
target[((hsi*150+143)*1+lsi)*1] = fp217;
LIBINT2_REALTYPE fp220;
fp220 = fp147 * fp138;
LIBINT2_REALTYPE fp219;
fp219 = fp220 * fp2;
target[((hsi*150+142)*1+lsi)*1] = fp219;
LIBINT2_REALTYPE fp242;
fp242 = fp147 * fp105;
LIBINT2_REALTYPE fp241;
fp241 = fp242 * fp7;
target[((hsi*150+131)*1+lsi)*1] = fp241;
LIBINT2_REALTYPE fp224;
fp224 = fp151 * fp138;
LIBINT2_REALTYPE fp223;
fp223 = fp224 * fp3;
target[((hsi*150+140)*1+lsi)*1] = fp223;
LIBINT2_REALTYPE fp4;
fp4 = 0 + src2[((hsi*35+24)*1+lsi)*1];
LIBINT2_REALTYPE fp226;
fp226 = fp141 * fp114;
LIBINT2_REALTYPE fp225;
fp225 = fp226 * fp4;
target[((hsi*150+139)*1+lsi)*1] = fp225;
LIBINT2_REALTYPE fp5;
fp5 = 0 + src2[((hsi*35+23)*1+lsi)*1];
LIBINT2_REALTYPE fp228;
fp228 = fp141 * fp105;
LIBINT2_REALTYPE fp227;
fp227 = fp228 * fp5;
target[((hsi*150+138)*1+lsi)*1] = fp227;
LIBINT2_REALTYPE fp6;
fp6 = 0 + src2[((hsi*35+22)*1+lsi)*1];
LIBINT2_REALTYPE fp230;
fp230 = fp141 * fp100;
LIBINT2_REALTYPE fp229;
fp229 = fp230 * fp6;
target[((hsi*150+137)*1+lsi)*1] = fp229;
LIBINT2_REALTYPE fp232;
fp232 = fp141 * fp95;
LIBINT2_REALTYPE fp231;
fp231 = fp232 * fp7;
target[((hsi*150+136)*1+lsi)*1] = fp231;
LIBINT2_REALTYPE fp234;
fp234 = fp144 * fp114;
LIBINT2_REALTYPE fp233;
fp233 = fp234 * fp5;
target[((hsi*150+135)*1+lsi)*1] = fp233;
LIBINT2_REALTYPE fp236;
fp236 = fp144 * fp105;
LIBINT2_REALTYPE fp235;
fp235 = fp236 * fp6;
target[((hsi*150+134)*1+lsi)*1] = fp235;
LIBINT2_REALTYPE fp238;
fp238 = fp144 * fp100;
LIBINT2_REALTYPE fp237;
fp237 = fp238 * fp7;
target[((hsi*150+133)*1+lsi)*1] = fp237;
LIBINT2_REALTYPE fp240;
fp240 = fp147 * fp114;
LIBINT2_REALTYPE fp239;
fp239 = fp240 * fp6;
target[((hsi*150+132)*1+lsi)*1] = fp239;
LIBINT2_REALTYPE fp222;
fp222 = fp147 * fp129;
LIBINT2_REALTYPE fp221;
fp221 = fp222 * fp3;
target[((hsi*150+141)*1+lsi)*1] = fp221;
LIBINT2_REALTYPE fp318;
fp318 = fp157 * fp124;
LIBINT2_REALTYPE fp317;
fp317 = fp318 * fp7;
target[((hsi*150+93)*1+lsi)*1] = fp317;
LIBINT2_REALTYPE fp320;
fp320 = fp160 * fp138;
LIBINT2_REALTYPE fp319;
fp319 = fp320 * fp6;
target[((hsi*150+92)*1+lsi)*1] = fp319;
LIBINT2_REALTYPE fp322;
fp322 = fp160 * fp129;
LIBINT2_REALTYPE fp321;
fp321 = fp322 * fp7;
target[((hsi*150+91)*1+lsi)*1] = fp321;
LIBINT2_REALTYPE fp324;
fp324 = fp164 * fp138;
LIBINT2_REALTYPE fp323;
fp323 = fp324 * fp7;
target[((hsi*150+90)*1+lsi)*1] = fp323;
LIBINT2_REALTYPE fp326;
fp326 = fp154 * fp114;
LIBINT2_REALTYPE fp325;
fp325 = fp326 * fp8;
target[((hsi*150+89)*1+lsi)*1] = fp325;
LIBINT2_REALTYPE fp328;
fp328 = fp154 * fp105;
LIBINT2_REALTYPE fp327;
fp327 = fp328 * fp9;
target[((hsi*150+88)*1+lsi)*1] = fp327;
LIBINT2_REALTYPE fp330;
fp330 = fp154 * fp100;
LIBINT2_REALTYPE fp329;
fp329 = fp330 * fp10;
target[((hsi*150+87)*1+lsi)*1] = fp329;
LIBINT2_REALTYPE fp332;
fp332 = fp154 * fp95;
LIBINT2_REALTYPE fp331;
fp331 = fp332 * fp11;
target[((hsi*150+86)*1+lsi)*1] = fp331;
LIBINT2_REALTYPE fp278;
fp278 = fp144 * fp52;
LIBINT2_REALTYPE fp277;
fp277 = fp278 * fp15;
target[((hsi*150+113)*1+lsi)*1] = fp277;
LIBINT2_REALTYPE fp336;
fp336 = fp157 * fp105;
LIBINT2_REALTYPE fp335;
fp335 = fp336 * fp10;
target[((hsi*150+84)*1+lsi)*1] = fp335;
LIBINT2_REALTYPE fp338;
fp338 = fp157 * fp100;
LIBINT2_REALTYPE fp337;
fp337 = fp338 * fp11;
target[((hsi*150+83)*1+lsi)*1] = fp337;
LIBINT2_REALTYPE fp340;
fp340 = fp160 * fp114;
LIBINT2_REALTYPE fp339;
fp339 = fp340 * fp10;
target[((hsi*150+82)*1+lsi)*1] = fp339;
LIBINT2_REALTYPE fp342;
fp342 = fp160 * fp105;
LIBINT2_REALTYPE fp341;
fp341 = fp342 * fp11;
target[((hsi*150+81)*1+lsi)*1] = fp341;
LIBINT2_REALTYPE fp344;
fp344 = fp164 * fp114;
LIBINT2_REALTYPE fp343;
fp343 = fp344 * fp11;
target[((hsi*150+80)*1+lsi)*1] = fp343;
LIBINT2_REALTYPE fp346;
fp346 = fp154 * fp90;
LIBINT2_REALTYPE fp345;
fp345 = fp346 * fp12;
target[((hsi*150+79)*1+lsi)*1] = fp345;
LIBINT2_REALTYPE fp348;
fp348 = fp154 * fp81;
LIBINT2_REALTYPE fp347;
fp347 = fp348 * fp13;
target[((hsi*150+78)*1+lsi)*1] = fp347;
LIBINT2_REALTYPE fp350;
fp350 = fp154 * fp76;
LIBINT2_REALTYPE fp349;
fp349 = fp350 * fp14;
target[((hsi*150+77)*1+lsi)*1] = fp349;
LIBINT2_REALTYPE fp352;
fp352 = fp154 * fp71;
LIBINT2_REALTYPE fp351;
fp351 = fp352 * fp15;
target[((hsi*150+76)*1+lsi)*1] = fp351;
LIBINT2_REALTYPE fp334;
fp334 = fp157 * fp114;
LIBINT2_REALTYPE fp333;
fp333 = fp334 * fp9;
target[((hsi*150+85)*1+lsi)*1] = fp333;
LIBINT2_REALTYPE fp280;
fp280 = fp147 * fp66;
LIBINT2_REALTYPE fp279;
fp279 = fp280 * fp14;
target[((hsi*150+112)*1+lsi)*1] = fp279;
LIBINT2_REALTYPE fp282;
fp282 = fp147 * fp57;
LIBINT2_REALTYPE fp281;
fp281 = fp282 * fp15;
target[((hsi*150+111)*1+lsi)*1] = fp281;
LIBINT2_REALTYPE fp284;
fp284 = fp151 * fp66;
LIBINT2_REALTYPE fp283;
fp283 = fp284 * fp15;
target[((hsi*150+110)*1+lsi)*1] = fp283;
LIBINT2_REALTYPE fp37;
fp37 = 0 + src1[((hsi*35+30)*1+lsi)*1];
LIBINT2_REALTYPE fp36;
fp36 = 0 + src1[((hsi*35+29)*1+lsi)*1];
LIBINT2_REALTYPE fp41;
fp41 = inteval->BO_y[vi] * fp36;
LIBINT2_REALTYPE fp40;
fp40 = fp37 + fp41;
LIBINT2_REALTYPE fp35;
fp35 = 0 + src1[((hsi*35+28)*1+lsi)*1];
LIBINT2_REALTYPE fp39;
fp39 = inteval->BO_y[vi] * fp35;
LIBINT2_REALTYPE fp38;
fp38 = fp36 + fp39;
LIBINT2_REALTYPE fp43;
fp43 = inteval->BO_y[vi] * fp38;
LIBINT2_REALTYPE fp42;
fp42 = fp40 + fp43;
LIBINT2_REALTYPE fp286;
fp286 = fp141 * fp42;
LIBINT2_REALTYPE fp285;
fp285 = fp286 * fp16;
target[((hsi*150+109)*1+lsi)*1] = fp285;
LIBINT2_REALTYPE fp30;
fp30 = 0 + src1[((hsi*35+31)*1+lsi)*1];
LIBINT2_REALTYPE fp32;
fp32 = inteval->BO_y[vi] * fp37;
LIBINT2_REALTYPE fp31;
fp31 = fp30 + fp32;
LIBINT2_REALTYPE fp34;
fp34 = inteval->BO_y[vi] * fp40;
LIBINT2_REALTYPE fp33;
fp33 = fp31 + fp34;
LIBINT2_REALTYPE fp288;
fp288 = fp141 * fp33;
LIBINT2_REALTYPE fp287;
fp287 = fp288 * fp17;
target[((hsi*150+108)*1+lsi)*1] = fp287;
LIBINT2_REALTYPE fp25;
fp25 = 0 + src1[((hsi*35+32)*1+lsi)*1];
LIBINT2_REALTYPE fp27;
fp27 = inteval->BO_y[vi] * fp30;
LIBINT2_REALTYPE fp26;
fp26 = fp25 + fp27;
LIBINT2_REALTYPE fp29;
fp29 = inteval->BO_y[vi] * fp31;
LIBINT2_REALTYPE fp28;
fp28 = fp26 + fp29;
LIBINT2_REALTYPE fp290;
fp290 = fp141 * fp28;
LIBINT2_REALTYPE fp289;
fp289 = fp290 * fp18;
target[((hsi*150+107)*1+lsi)*1] = fp289;
LIBINT2_REALTYPE fp24;
fp24 = inteval->BO_y[vi] * fp26;
LIBINT2_REALTYPE fp22;
fp22 = inteval->BO_y[vi] * fp25;
LIBINT2_REALTYPE fp20;
fp20 = 0 + src1[((hsi*35+33)*1+lsi)*1];
LIBINT2_REALTYPE fp21;
fp21 = fp20 + fp22;
LIBINT2_REALTYPE fp23;
fp23 = fp21 + fp24;
LIBINT2_REALTYPE fp292;
fp292 = fp141 * fp23;
LIBINT2_REALTYPE fp291;
fp291 = fp292 * fp19;
target[((hsi*150+106)*1+lsi)*1] = fp291;
LIBINT2_REALTYPE fp294;
fp294 = fp144 * fp42;
LIBINT2_REALTYPE fp293;
fp293 = fp294 * fp17;
target[((hsi*150+105)*1+lsi)*1] = fp293;
LIBINT2_REALTYPE fp316;
fp316 = fp157 * fp129;
LIBINT2_REALTYPE fp315;
fp315 = fp316 * fp6;
target[((hsi*150+94)*1+lsi)*1] = fp315;
LIBINT2_REALTYPE fp298;
fp298 = fp144 * fp28;
LIBINT2_REALTYPE fp297;
fp297 = fp298 * fp19;
target[((hsi*150+103)*1+lsi)*1] = fp297;
LIBINT2_REALTYPE fp300;
fp300 = fp147 * fp42;
LIBINT2_REALTYPE fp299;
fp299 = fp300 * fp18;
target[((hsi*150+102)*1+lsi)*1] = fp299;
LIBINT2_REALTYPE fp302;
fp302 = fp147 * fp33;
LIBINT2_REALTYPE fp301;
fp301 = fp302 * fp19;
target[((hsi*150+101)*1+lsi)*1] = fp301;
LIBINT2_REALTYPE fp304;
fp304 = fp151 * fp42;
LIBINT2_REALTYPE fp303;
fp303 = fp304 * fp19;
target[((hsi*150+100)*1+lsi)*1] = fp303;
LIBINT2_REALTYPE fp306;
fp306 = fp154 * fp138;
LIBINT2_REALTYPE fp305;
fp305 = fp306 * fp4;
target[((hsi*150+99)*1+lsi)*1] = fp305;
LIBINT2_REALTYPE fp308;
fp308 = fp154 * fp129;
LIBINT2_REALTYPE fp307;
fp307 = fp308 * fp5;
target[((hsi*150+98)*1+lsi)*1] = fp307;
LIBINT2_REALTYPE fp310;
fp310 = fp154 * fp124;
LIBINT2_REALTYPE fp309;
fp309 = fp310 * fp6;
target[((hsi*150+97)*1+lsi)*1] = fp309;
LIBINT2_REALTYPE fp312;
fp312 = fp154 * fp119;
LIBINT2_REALTYPE fp311;
fp311 = fp312 * fp7;
target[((hsi*150+96)*1+lsi)*1] = fp311;
LIBINT2_REALTYPE fp314;
fp314 = fp157 * fp138;
LIBINT2_REALTYPE fp313;
fp313 = fp314 * fp5;
target[((hsi*150+95)*1+lsi)*1] = fp313;
LIBINT2_REALTYPE fp296;
fp296 = fp144 * fp33;
LIBINT2_REALTYPE fp295;
fp295 = fp296 * fp18;
target[((hsi*150+104)*1+lsi)*1] = fp295;
}
}
}
const int hsi = 0;
const int lsi = 0;
const int vi = 0;
/** Number of flops = 505 */
}

#ifdef __cplusplus
};
#endif
