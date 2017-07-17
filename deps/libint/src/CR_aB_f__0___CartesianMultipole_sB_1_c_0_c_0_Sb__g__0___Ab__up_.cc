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
void CR_aB_f__0___CartesianMultipole_sB_1_c_0_c_0_Sb__g__0___Ab__up_(const Libint_t* inteval, LIBINT2_REALTYPE* target, const LIBINT2_REALTYPE* src0, const LIBINT2_REALTYPE* src1, const LIBINT2_REALTYPE* src2) {

LIBINT2_REALTYPE*const  stack = target;
{
const int hsi = 0;
{
const int lsi = 0;
{
const int vi = 0;
LIBINT2_REALTYPE fp19;
fp19 = 0 + src2[((hsi*24+0)*1+lsi)*1];
LIBINT2_REALTYPE fp59;
fp59 = 0 + src0[((hsi*24+7)*1+lsi)*1];
LIBINT2_REALTYPE fp56;
fp56 = 0 + src0[((hsi*24+6)*1+lsi)*1];
LIBINT2_REALTYPE fp58;
fp58 = inteval->BO_x[vi] * fp56;
LIBINT2_REALTYPE fp57;
fp57 = fp59 + fp58;
LIBINT2_REALTYPE fp25;
fp25 = 0 + src1[((hsi*24+16)*1+lsi)*1];
LIBINT2_REALTYPE fp293;
fp293 = fp57 * fp25;
LIBINT2_REALTYPE fp292;
fp292 = fp293 * fp19;
target[((hsi*150+55)*1+lsi)*1] = fp292;
LIBINT2_REALTYPE fp16;
fp16 = 0 + src2[((hsi*24+3)*1+lsi)*1];
LIBINT2_REALTYPE fp62;
fp62 = 0 + src0[((hsi*24+8)*1+lsi)*1];
LIBINT2_REALTYPE fp61;
fp61 = inteval->BO_x[vi] * fp59;
LIBINT2_REALTYPE fp60;
fp60 = fp62 + fp61;
LIBINT2_REALTYPE fp29;
fp29 = 0 + src1[((hsi*24+12)*1+lsi)*1];
LIBINT2_REALTYPE fp295;
fp295 = fp60 * fp29;
LIBINT2_REALTYPE fp294;
fp294 = fp295 * fp16;
target[((hsi*150+54)*1+lsi)*1] = fp294;
LIBINT2_REALTYPE fp17;
fp17 = 0 + src2[((hsi*24+2)*1+lsi)*1];
LIBINT2_REALTYPE fp28;
fp28 = 0 + src1[((hsi*24+13)*1+lsi)*1];
LIBINT2_REALTYPE fp297;
fp297 = fp60 * fp28;
LIBINT2_REALTYPE fp296;
fp296 = fp297 * fp17;
target[((hsi*150+53)*1+lsi)*1] = fp296;
LIBINT2_REALTYPE fp18;
fp18 = 0 + src2[((hsi*24+1)*1+lsi)*1];
LIBINT2_REALTYPE fp27;
fp27 = 0 + src1[((hsi*24+14)*1+lsi)*1];
LIBINT2_REALTYPE fp299;
fp299 = fp60 * fp27;
LIBINT2_REALTYPE fp298;
fp298 = fp299 * fp18;
target[((hsi*150+52)*1+lsi)*1] = fp298;
LIBINT2_REALTYPE fp26;
fp26 = 0 + src1[((hsi*24+15)*1+lsi)*1];
LIBINT2_REALTYPE fp301;
fp301 = fp60 * fp26;
LIBINT2_REALTYPE fp300;
fp300 = fp301 * fp19;
target[((hsi*150+51)*1+lsi)*1] = fp300;
LIBINT2_REALTYPE fp65;
fp65 = 0 + src0[((hsi*24+9)*1+lsi)*1];
LIBINT2_REALTYPE fp64;
fp64 = inteval->BO_x[vi] * fp62;
LIBINT2_REALTYPE fp63;
fp63 = fp65 + fp64;
LIBINT2_REALTYPE fp303;
fp303 = fp63 * fp29;
LIBINT2_REALTYPE fp302;
fp302 = fp303 * fp17;
target[((hsi*150+50)*1+lsi)*1] = fp302;
LIBINT2_REALTYPE fp305;
fp305 = fp63 * fp28;
LIBINT2_REALTYPE fp304;
fp304 = fp305 * fp18;
target[((hsi*150+49)*1+lsi)*1] = fp304;
LIBINT2_REALTYPE fp307;
fp307 = fp63 * fp27;
LIBINT2_REALTYPE fp306;
fp306 = fp307 * fp19;
target[((hsi*150+48)*1+lsi)*1] = fp306;
LIBINT2_REALTYPE fp39;
fp39 = 0 + src1[((hsi*24+0)*1+lsi)*1];
LIBINT2_REALTYPE fp100;
fp100 = 0 + src0[((hsi*24+22)*1+lsi)*1];
LIBINT2_REALTYPE fp103;
fp103 = inteval->BO_x[vi] * fp100;
LIBINT2_REALTYPE fp101;
fp101 = 0 + src0[((hsi*24+23)*1+lsi)*1];
LIBINT2_REALTYPE fp102;
fp102 = fp101 + fp103;
LIBINT2_REALTYPE fp403;
fp403 = fp102 * fp39;
LIBINT2_REALTYPE fp402;
fp402 = fp403 * fp19;
target[((hsi*150+0)*1+lsi)*1] = fp402;
LIBINT2_REALTYPE fp68;
fp68 = 0 + src0[((hsi*24+10)*1+lsi)*1];
LIBINT2_REALTYPE fp67;
fp67 = inteval->BO_x[vi] * fp65;
LIBINT2_REALTYPE fp66;
fp66 = fp68 + fp67;
LIBINT2_REALTYPE fp311;
fp311 = fp66 * fp28;
LIBINT2_REALTYPE fp310;
fp310 = fp311 * fp19;
target[((hsi*150+46)*1+lsi)*1] = fp310;
LIBINT2_REALTYPE fp71;
fp71 = inteval->BO_x[vi] * fp68;
LIBINT2_REALTYPE fp69;
fp69 = 0 + src0[((hsi*24+11)*1+lsi)*1];
LIBINT2_REALTYPE fp70;
fp70 = fp69 + fp71;
LIBINT2_REALTYPE fp313;
fp313 = fp70 * fp29;
LIBINT2_REALTYPE fp312;
fp312 = fp313 * fp19;
target[((hsi*150+45)*1+lsi)*1] = fp312;
LIBINT2_REALTYPE fp10;
fp10 = 0 + src2[((hsi*24+10)*1+lsi)*1];
LIBINT2_REALTYPE fp75;
fp75 = 0 + src0[((hsi*24+13)*1+lsi)*1];
LIBINT2_REALTYPE fp72;
fp72 = 0 + src0[((hsi*24+12)*1+lsi)*1];
LIBINT2_REALTYPE fp74;
fp74 = inteval->BO_x[vi] * fp72;
LIBINT2_REALTYPE fp73;
fp73 = fp75 + fp74;
LIBINT2_REALTYPE fp315;
fp315 = fp73 * fp39;
LIBINT2_REALTYPE fp314;
fp314 = fp315 * fp10;
target[((hsi*150+44)*1+lsi)*1] = fp314;
LIBINT2_REALTYPE fp11;
fp11 = 0 + src2[((hsi*24+9)*1+lsi)*1];
LIBINT2_REALTYPE fp38;
fp38 = 0 + src1[((hsi*24+1)*1+lsi)*1];
LIBINT2_REALTYPE fp317;
fp317 = fp73 * fp38;
LIBINT2_REALTYPE fp316;
fp316 = fp317 * fp11;
target[((hsi*150+43)*1+lsi)*1] = fp316;
LIBINT2_REALTYPE fp12;
fp12 = 0 + src2[((hsi*24+8)*1+lsi)*1];
LIBINT2_REALTYPE fp37;
fp37 = 0 + src1[((hsi*24+2)*1+lsi)*1];
LIBINT2_REALTYPE fp319;
fp319 = fp73 * fp37;
LIBINT2_REALTYPE fp318;
fp318 = fp319 * fp12;
target[((hsi*150+42)*1+lsi)*1] = fp318;
LIBINT2_REALTYPE fp13;
fp13 = 0 + src2[((hsi*24+7)*1+lsi)*1];
LIBINT2_REALTYPE fp36;
fp36 = 0 + src1[((hsi*24+3)*1+lsi)*1];
LIBINT2_REALTYPE fp321;
fp321 = fp73 * fp36;
LIBINT2_REALTYPE fp320;
fp320 = fp321 * fp13;
target[((hsi*150+41)*1+lsi)*1] = fp320;
LIBINT2_REALTYPE fp14;
fp14 = 0 + src2[((hsi*24+6)*1+lsi)*1];
LIBINT2_REALTYPE fp35;
fp35 = 0 + src1[((hsi*24+4)*1+lsi)*1];
LIBINT2_REALTYPE fp323;
fp323 = fp73 * fp35;
LIBINT2_REALTYPE fp322;
fp322 = fp323 * fp14;
target[((hsi*150+40)*1+lsi)*1] = fp322;
LIBINT2_REALTYPE fp78;
fp78 = 0 + src0[((hsi*24+14)*1+lsi)*1];
LIBINT2_REALTYPE fp77;
fp77 = inteval->BO_x[vi] * fp75;
LIBINT2_REALTYPE fp76;
fp76 = fp78 + fp77;
LIBINT2_REALTYPE fp325;
fp325 = fp76 * fp39;
LIBINT2_REALTYPE fp324;
fp324 = fp325 * fp11;
target[((hsi*150+39)*1+lsi)*1] = fp324;
LIBINT2_REALTYPE fp309;
fp309 = fp66 * fp29;
LIBINT2_REALTYPE fp308;
fp308 = fp309 * fp18;
target[((hsi*150+47)*1+lsi)*1] = fp308;
LIBINT2_REALTYPE fp34;
fp34 = 0 + src1[((hsi*24+6)*1+lsi)*1];
LIBINT2_REALTYPE fp255;
fp255 = fp57 * fp34;
LIBINT2_REALTYPE fp254;
fp254 = fp255 * fp10;
target[((hsi*150+74)*1+lsi)*1] = fp254;
LIBINT2_REALTYPE fp33;
fp33 = 0 + src1[((hsi*24+7)*1+lsi)*1];
LIBINT2_REALTYPE fp257;
fp257 = fp57 * fp33;
LIBINT2_REALTYPE fp256;
fp256 = fp257 * fp11;
target[((hsi*150+73)*1+lsi)*1] = fp256;
LIBINT2_REALTYPE fp32;
fp32 = 0 + src1[((hsi*24+8)*1+lsi)*1];
LIBINT2_REALTYPE fp259;
fp259 = fp57 * fp32;
LIBINT2_REALTYPE fp258;
fp258 = fp259 * fp12;
target[((hsi*150+72)*1+lsi)*1] = fp258;
LIBINT2_REALTYPE fp31;
fp31 = 0 + src1[((hsi*24+9)*1+lsi)*1];
LIBINT2_REALTYPE fp261;
fp261 = fp57 * fp31;
LIBINT2_REALTYPE fp260;
fp260 = fp261 * fp13;
target[((hsi*150+71)*1+lsi)*1] = fp260;
LIBINT2_REALTYPE fp30;
fp30 = 0 + src1[((hsi*24+10)*1+lsi)*1];
LIBINT2_REALTYPE fp263;
fp263 = fp57 * fp30;
LIBINT2_REALTYPE fp262;
fp262 = fp263 * fp14;
target[((hsi*150+70)*1+lsi)*1] = fp262;
LIBINT2_REALTYPE fp265;
fp265 = fp60 * fp34;
LIBINT2_REALTYPE fp264;
fp264 = fp265 * fp11;
target[((hsi*150+69)*1+lsi)*1] = fp264;
LIBINT2_REALTYPE fp267;
fp267 = fp60 * fp33;
LIBINT2_REALTYPE fp266;
fp266 = fp267 * fp12;
target[((hsi*150+68)*1+lsi)*1] = fp266;
LIBINT2_REALTYPE fp269;
fp269 = fp60 * fp32;
LIBINT2_REALTYPE fp268;
fp268 = fp269 * fp13;
target[((hsi*150+67)*1+lsi)*1] = fp268;
LIBINT2_REALTYPE fp291;
fp291 = fp57 * fp26;
LIBINT2_REALTYPE fp290;
fp290 = fp291 * fp18;
target[((hsi*150+56)*1+lsi)*1] = fp290;
LIBINT2_REALTYPE fp273;
fp273 = fp63 * fp34;
LIBINT2_REALTYPE fp272;
fp272 = fp273 * fp12;
target[((hsi*150+65)*1+lsi)*1] = fp272;
LIBINT2_REALTYPE fp275;
fp275 = fp63 * fp33;
LIBINT2_REALTYPE fp274;
fp274 = fp275 * fp13;
target[((hsi*150+64)*1+lsi)*1] = fp274;
LIBINT2_REALTYPE fp277;
fp277 = fp63 * fp32;
LIBINT2_REALTYPE fp276;
fp276 = fp277 * fp14;
target[((hsi*150+63)*1+lsi)*1] = fp276;
LIBINT2_REALTYPE fp279;
fp279 = fp66 * fp34;
LIBINT2_REALTYPE fp278;
fp278 = fp279 * fp13;
target[((hsi*150+62)*1+lsi)*1] = fp278;
LIBINT2_REALTYPE fp281;
fp281 = fp66 * fp33;
LIBINT2_REALTYPE fp280;
fp280 = fp281 * fp14;
target[((hsi*150+61)*1+lsi)*1] = fp280;
LIBINT2_REALTYPE fp283;
fp283 = fp70 * fp34;
LIBINT2_REALTYPE fp282;
fp282 = fp283 * fp14;
target[((hsi*150+60)*1+lsi)*1] = fp282;
LIBINT2_REALTYPE fp15;
fp15 = 0 + src2[((hsi*24+4)*1+lsi)*1];
LIBINT2_REALTYPE fp285;
fp285 = fp57 * fp29;
LIBINT2_REALTYPE fp284;
fp284 = fp285 * fp15;
target[((hsi*150+59)*1+lsi)*1] = fp284;
LIBINT2_REALTYPE fp287;
fp287 = fp57 * fp28;
LIBINT2_REALTYPE fp286;
fp286 = fp287 * fp16;
target[((hsi*150+58)*1+lsi)*1] = fp286;
LIBINT2_REALTYPE fp289;
fp289 = fp57 * fp27;
LIBINT2_REALTYPE fp288;
fp288 = fp289 * fp17;
target[((hsi*150+57)*1+lsi)*1] = fp288;
LIBINT2_REALTYPE fp271;
fp271 = fp60 * fp31;
LIBINT2_REALTYPE fp270;
fp270 = fp271 * fp14;
target[((hsi*150+66)*1+lsi)*1] = fp270;
LIBINT2_REALTYPE fp81;
fp81 = 0 + src0[((hsi*24+15)*1+lsi)*1];
LIBINT2_REALTYPE fp80;
fp80 = inteval->BO_x[vi] * fp78;
LIBINT2_REALTYPE fp79;
fp79 = fp81 + fp80;
LIBINT2_REALTYPE fp367;
fp367 = fp79 * fp32;
LIBINT2_REALTYPE fp366;
fp366 = fp367 * fp19;
target[((hsi*150+18)*1+lsi)*1] = fp366;
LIBINT2_REALTYPE fp84;
fp84 = 0 + src0[((hsi*24+16)*1+lsi)*1];
LIBINT2_REALTYPE fp83;
fp83 = inteval->BO_x[vi] * fp81;
LIBINT2_REALTYPE fp82;
fp82 = fp84 + fp83;
LIBINT2_REALTYPE fp369;
fp369 = fp82 * fp34;
LIBINT2_REALTYPE fp368;
fp368 = fp369 * fp18;
target[((hsi*150+17)*1+lsi)*1] = fp368;
LIBINT2_REALTYPE fp371;
fp371 = fp82 * fp33;
LIBINT2_REALTYPE fp370;
fp370 = fp371 * fp19;
target[((hsi*150+16)*1+lsi)*1] = fp370;
LIBINT2_REALTYPE fp87;
fp87 = inteval->BO_x[vi] * fp84;
LIBINT2_REALTYPE fp85;
fp85 = 0 + src0[((hsi*24+17)*1+lsi)*1];
LIBINT2_REALTYPE fp86;
fp86 = fp85 + fp87;
LIBINT2_REALTYPE fp373;
fp373 = fp86 * fp34;
LIBINT2_REALTYPE fp372;
fp372 = fp373 * fp19;
target[((hsi*150+15)*1+lsi)*1] = fp372;
LIBINT2_REALTYPE fp91;
fp91 = 0 + src0[((hsi*24+19)*1+lsi)*1];
LIBINT2_REALTYPE fp88;
fp88 = 0 + src0[((hsi*24+18)*1+lsi)*1];
LIBINT2_REALTYPE fp90;
fp90 = inteval->BO_x[vi] * fp88;
LIBINT2_REALTYPE fp89;
fp89 = fp91 + fp90;
LIBINT2_REALTYPE fp375;
fp375 = fp89 * fp39;
LIBINT2_REALTYPE fp374;
fp374 = fp375 * fp15;
target[((hsi*150+14)*1+lsi)*1] = fp374;
LIBINT2_REALTYPE fp377;
fp377 = fp89 * fp38;
LIBINT2_REALTYPE fp376;
fp376 = fp377 * fp16;
target[((hsi*150+13)*1+lsi)*1] = fp376;
LIBINT2_REALTYPE fp379;
fp379 = fp89 * fp37;
LIBINT2_REALTYPE fp378;
fp378 = fp379 * fp17;
target[((hsi*150+12)*1+lsi)*1] = fp378;
LIBINT2_REALTYPE fp381;
fp381 = fp89 * fp36;
LIBINT2_REALTYPE fp380;
fp380 = fp381 * fp18;
target[((hsi*150+11)*1+lsi)*1] = fp380;
LIBINT2_REALTYPE fp327;
fp327 = fp76 * fp38;
LIBINT2_REALTYPE fp326;
fp326 = fp327 * fp12;
target[((hsi*150+38)*1+lsi)*1] = fp326;
LIBINT2_REALTYPE fp94;
fp94 = 0 + src0[((hsi*24+20)*1+lsi)*1];
LIBINT2_REALTYPE fp93;
fp93 = inteval->BO_x[vi] * fp91;
LIBINT2_REALTYPE fp92;
fp92 = fp94 + fp93;
LIBINT2_REALTYPE fp385;
fp385 = fp92 * fp39;
LIBINT2_REALTYPE fp384;
fp384 = fp385 * fp16;
target[((hsi*150+9)*1+lsi)*1] = fp384;
LIBINT2_REALTYPE fp387;
fp387 = fp92 * fp38;
LIBINT2_REALTYPE fp386;
fp386 = fp387 * fp17;
target[((hsi*150+8)*1+lsi)*1] = fp386;
LIBINT2_REALTYPE fp389;
fp389 = fp92 * fp37;
LIBINT2_REALTYPE fp388;
fp388 = fp389 * fp18;
target[((hsi*150+7)*1+lsi)*1] = fp388;
LIBINT2_REALTYPE fp391;
fp391 = fp92 * fp36;
LIBINT2_REALTYPE fp390;
fp390 = fp391 * fp19;
target[((hsi*150+6)*1+lsi)*1] = fp390;
LIBINT2_REALTYPE fp97;
fp97 = 0 + src0[((hsi*24+21)*1+lsi)*1];
LIBINT2_REALTYPE fp96;
fp96 = inteval->BO_x[vi] * fp94;
LIBINT2_REALTYPE fp95;
fp95 = fp97 + fp96;
LIBINT2_REALTYPE fp393;
fp393 = fp95 * fp39;
LIBINT2_REALTYPE fp392;
fp392 = fp393 * fp17;
target[((hsi*150+5)*1+lsi)*1] = fp392;
LIBINT2_REALTYPE fp395;
fp395 = fp95 * fp38;
LIBINT2_REALTYPE fp394;
fp394 = fp395 * fp18;
target[((hsi*150+4)*1+lsi)*1] = fp394;
LIBINT2_REALTYPE fp397;
fp397 = fp95 * fp37;
LIBINT2_REALTYPE fp396;
fp396 = fp397 * fp19;
target[((hsi*150+3)*1+lsi)*1] = fp396;
LIBINT2_REALTYPE fp99;
fp99 = inteval->BO_x[vi] * fp97;
LIBINT2_REALTYPE fp98;
fp98 = fp100 + fp99;
LIBINT2_REALTYPE fp399;
fp399 = fp98 * fp39;
LIBINT2_REALTYPE fp398;
fp398 = fp399 * fp18;
target[((hsi*150+2)*1+lsi)*1] = fp398;
LIBINT2_REALTYPE fp401;
fp401 = fp98 * fp38;
LIBINT2_REALTYPE fp400;
fp400 = fp401 * fp19;
target[((hsi*150+1)*1+lsi)*1] = fp400;
LIBINT2_REALTYPE fp383;
fp383 = fp89 * fp35;
LIBINT2_REALTYPE fp382;
fp382 = fp383 * fp19;
target[((hsi*150+10)*1+lsi)*1] = fp382;
LIBINT2_REALTYPE fp329;
fp329 = fp76 * fp37;
LIBINT2_REALTYPE fp328;
fp328 = fp329 * fp13;
target[((hsi*150+37)*1+lsi)*1] = fp328;
LIBINT2_REALTYPE fp331;
fp331 = fp76 * fp36;
LIBINT2_REALTYPE fp330;
fp330 = fp331 * fp14;
target[((hsi*150+36)*1+lsi)*1] = fp330;
LIBINT2_REALTYPE fp333;
fp333 = fp79 * fp39;
LIBINT2_REALTYPE fp332;
fp332 = fp333 * fp12;
target[((hsi*150+35)*1+lsi)*1] = fp332;
LIBINT2_REALTYPE fp335;
fp335 = fp79 * fp38;
LIBINT2_REALTYPE fp334;
fp334 = fp335 * fp13;
target[((hsi*150+34)*1+lsi)*1] = fp334;
LIBINT2_REALTYPE fp337;
fp337 = fp79 * fp37;
LIBINT2_REALTYPE fp336;
fp336 = fp337 * fp14;
target[((hsi*150+33)*1+lsi)*1] = fp336;
LIBINT2_REALTYPE fp339;
fp339 = fp82 * fp39;
LIBINT2_REALTYPE fp338;
fp338 = fp339 * fp13;
target[((hsi*150+32)*1+lsi)*1] = fp338;
LIBINT2_REALTYPE fp341;
fp341 = fp82 * fp38;
LIBINT2_REALTYPE fp340;
fp340 = fp341 * fp14;
target[((hsi*150+31)*1+lsi)*1] = fp340;
LIBINT2_REALTYPE fp343;
fp343 = fp86 * fp39;
LIBINT2_REALTYPE fp342;
fp342 = fp343 * fp14;
target[((hsi*150+30)*1+lsi)*1] = fp342;
LIBINT2_REALTYPE fp365;
fp365 = fp79 * fp33;
LIBINT2_REALTYPE fp364;
fp364 = fp365 * fp18;
target[((hsi*150+19)*1+lsi)*1] = fp364;
LIBINT2_REALTYPE fp347;
fp347 = fp73 * fp33;
LIBINT2_REALTYPE fp346;
fp346 = fp347 * fp16;
target[((hsi*150+28)*1+lsi)*1] = fp346;
LIBINT2_REALTYPE fp349;
fp349 = fp73 * fp32;
LIBINT2_REALTYPE fp348;
fp348 = fp349 * fp17;
target[((hsi*150+27)*1+lsi)*1] = fp348;
LIBINT2_REALTYPE fp351;
fp351 = fp73 * fp31;
LIBINT2_REALTYPE fp350;
fp350 = fp351 * fp18;
target[((hsi*150+26)*1+lsi)*1] = fp350;
LIBINT2_REALTYPE fp353;
fp353 = fp73 * fp30;
LIBINT2_REALTYPE fp352;
fp352 = fp353 * fp19;
target[((hsi*150+25)*1+lsi)*1] = fp352;
LIBINT2_REALTYPE fp355;
fp355 = fp76 * fp34;
LIBINT2_REALTYPE fp354;
fp354 = fp355 * fp16;
target[((hsi*150+24)*1+lsi)*1] = fp354;
LIBINT2_REALTYPE fp357;
fp357 = fp76 * fp33;
LIBINT2_REALTYPE fp356;
fp356 = fp357 * fp17;
target[((hsi*150+23)*1+lsi)*1] = fp356;
LIBINT2_REALTYPE fp359;
fp359 = fp76 * fp32;
LIBINT2_REALTYPE fp358;
fp358 = fp359 * fp18;
target[((hsi*150+22)*1+lsi)*1] = fp358;
LIBINT2_REALTYPE fp361;
fp361 = fp76 * fp31;
LIBINT2_REALTYPE fp360;
fp360 = fp361 * fp19;
target[((hsi*150+21)*1+lsi)*1] = fp360;
LIBINT2_REALTYPE fp363;
fp363 = fp79 * fp34;
LIBINT2_REALTYPE fp362;
fp362 = fp363 * fp17;
target[((hsi*150+20)*1+lsi)*1] = fp362;
LIBINT2_REALTYPE fp345;
fp345 = fp73 * fp34;
LIBINT2_REALTYPE fp344;
fp344 = fp345 * fp15;
target[((hsi*150+29)*1+lsi)*1] = fp344;
LIBINT2_REALTYPE fp9;
fp9 = 0 + src2[((hsi*24+12)*1+lsi)*1];
LIBINT2_REALTYPE fp43;
fp43 = 0 + src0[((hsi*24+1)*1+lsi)*1];
LIBINT2_REALTYPE fp40;
fp40 = 0 + src0[((hsi*24+0)*1+lsi)*1];
LIBINT2_REALTYPE fp42;
fp42 = inteval->BO_x[vi] * fp40;
LIBINT2_REALTYPE fp41;
fp41 = fp43 + fp42;
LIBINT2_REALTYPE fp143;
fp143 = fp41 * fp30;
LIBINT2_REALTYPE fp142;
fp142 = fp143 * fp9;
target[((hsi*150+130)*1+lsi)*1] = fp142;
LIBINT2_REALTYPE fp6;
fp6 = 0 + src2[((hsi*24+15)*1+lsi)*1];
LIBINT2_REALTYPE fp46;
fp46 = 0 + src0[((hsi*24+2)*1+lsi)*1];
LIBINT2_REALTYPE fp45;
fp45 = inteval->BO_x[vi] * fp43;
LIBINT2_REALTYPE fp44;
fp44 = fp46 + fp45;
LIBINT2_REALTYPE fp145;
fp145 = fp44 * fp34;
LIBINT2_REALTYPE fp144;
fp144 = fp145 * fp6;
target[((hsi*150+129)*1+lsi)*1] = fp144;
LIBINT2_REALTYPE fp7;
fp7 = 0 + src2[((hsi*24+14)*1+lsi)*1];
LIBINT2_REALTYPE fp147;
fp147 = fp44 * fp33;
LIBINT2_REALTYPE fp146;
fp146 = fp147 * fp7;
target[((hsi*150+128)*1+lsi)*1] = fp146;
LIBINT2_REALTYPE fp8;
fp8 = 0 + src2[((hsi*24+13)*1+lsi)*1];
LIBINT2_REALTYPE fp149;
fp149 = fp44 * fp32;
LIBINT2_REALTYPE fp148;
fp148 = fp149 * fp8;
target[((hsi*150+127)*1+lsi)*1] = fp148;
LIBINT2_REALTYPE fp151;
fp151 = fp44 * fp31;
LIBINT2_REALTYPE fp150;
fp150 = fp151 * fp9;
target[((hsi*150+126)*1+lsi)*1] = fp150;
LIBINT2_REALTYPE fp49;
fp49 = 0 + src0[((hsi*24+3)*1+lsi)*1];
LIBINT2_REALTYPE fp48;
fp48 = inteval->BO_x[vi] * fp46;
LIBINT2_REALTYPE fp47;
fp47 = fp49 + fp48;
LIBINT2_REALTYPE fp153;
fp153 = fp47 * fp34;
LIBINT2_REALTYPE fp152;
fp152 = fp153 * fp7;
target[((hsi*150+125)*1+lsi)*1] = fp152;
LIBINT2_REALTYPE fp155;
fp155 = fp47 * fp33;
LIBINT2_REALTYPE fp154;
fp154 = fp155 * fp8;
target[((hsi*150+124)*1+lsi)*1] = fp154;
LIBINT2_REALTYPE fp157;
fp157 = fp47 * fp32;
LIBINT2_REALTYPE fp156;
fp156 = fp157 * fp9;
target[((hsi*150+123)*1+lsi)*1] = fp156;
LIBINT2_REALTYPE fp253;
fp253 = fp70 * fp39;
LIBINT2_REALTYPE fp252;
fp252 = fp253 * fp9;
target[((hsi*150+75)*1+lsi)*1] = fp252;
LIBINT2_REALTYPE fp52;
fp52 = 0 + src0[((hsi*24+4)*1+lsi)*1];
LIBINT2_REALTYPE fp51;
fp51 = inteval->BO_x[vi] * fp49;
LIBINT2_REALTYPE fp50;
fp50 = fp52 + fp51;
LIBINT2_REALTYPE fp161;
fp161 = fp50 * fp33;
LIBINT2_REALTYPE fp160;
fp160 = fp161 * fp9;
target[((hsi*150+121)*1+lsi)*1] = fp160;
LIBINT2_REALTYPE fp55;
fp55 = inteval->BO_x[vi] * fp52;
LIBINT2_REALTYPE fp53;
fp53 = 0 + src0[((hsi*24+5)*1+lsi)*1];
LIBINT2_REALTYPE fp54;
fp54 = fp53 + fp55;
LIBINT2_REALTYPE fp163;
fp163 = fp54 * fp34;
LIBINT2_REALTYPE fp162;
fp162 = fp163 * fp9;
target[((hsi*150+120)*1+lsi)*1] = fp162;
LIBINT2_REALTYPE fp165;
fp165 = fp41 * fp29;
LIBINT2_REALTYPE fp164;
fp164 = fp165 * fp10;
target[((hsi*150+119)*1+lsi)*1] = fp164;
LIBINT2_REALTYPE fp167;
fp167 = fp41 * fp28;
LIBINT2_REALTYPE fp166;
fp166 = fp167 * fp11;
target[((hsi*150+118)*1+lsi)*1] = fp166;
LIBINT2_REALTYPE fp169;
fp169 = fp41 * fp27;
LIBINT2_REALTYPE fp168;
fp168 = fp169 * fp12;
target[((hsi*150+117)*1+lsi)*1] = fp168;
LIBINT2_REALTYPE fp171;
fp171 = fp41 * fp26;
LIBINT2_REALTYPE fp170;
fp170 = fp171 * fp13;
target[((hsi*150+116)*1+lsi)*1] = fp170;
LIBINT2_REALTYPE fp173;
fp173 = fp41 * fp25;
LIBINT2_REALTYPE fp172;
fp172 = fp173 * fp14;
target[((hsi*150+115)*1+lsi)*1] = fp172;
LIBINT2_REALTYPE fp175;
fp175 = fp44 * fp29;
LIBINT2_REALTYPE fp174;
fp174 = fp175 * fp11;
target[((hsi*150+114)*1+lsi)*1] = fp174;
LIBINT2_REALTYPE fp159;
fp159 = fp50 * fp34;
LIBINT2_REALTYPE fp158;
fp158 = fp159 * fp8;
target[((hsi*150+122)*1+lsi)*1] = fp158;
LIBINT2_REALTYPE fp0;
fp0 = 0 + src2[((hsi*24+22)*1+lsi)*1];
LIBINT2_REALTYPE fp105;
fp105 = fp41 * fp39;
LIBINT2_REALTYPE fp104;
fp104 = fp105 * fp0;
target[((hsi*150+149)*1+lsi)*1] = fp104;
LIBINT2_REALTYPE fp1;
fp1 = 0 + src2[((hsi*24+21)*1+lsi)*1];
LIBINT2_REALTYPE fp107;
fp107 = fp41 * fp38;
LIBINT2_REALTYPE fp106;
fp106 = fp107 * fp1;
target[((hsi*150+148)*1+lsi)*1] = fp106;
LIBINT2_REALTYPE fp2;
fp2 = 0 + src2[((hsi*24+20)*1+lsi)*1];
LIBINT2_REALTYPE fp109;
fp109 = fp41 * fp37;
LIBINT2_REALTYPE fp108;
fp108 = fp109 * fp2;
target[((hsi*150+147)*1+lsi)*1] = fp108;
LIBINT2_REALTYPE fp3;
fp3 = 0 + src2[((hsi*24+19)*1+lsi)*1];
LIBINT2_REALTYPE fp111;
fp111 = fp41 * fp36;
LIBINT2_REALTYPE fp110;
fp110 = fp111 * fp3;
target[((hsi*150+146)*1+lsi)*1] = fp110;
LIBINT2_REALTYPE fp4;
fp4 = 0 + src2[((hsi*24+18)*1+lsi)*1];
LIBINT2_REALTYPE fp113;
fp113 = fp41 * fp35;
LIBINT2_REALTYPE fp112;
fp112 = fp113 * fp4;
target[((hsi*150+145)*1+lsi)*1] = fp112;
LIBINT2_REALTYPE fp115;
fp115 = fp44 * fp39;
LIBINT2_REALTYPE fp114;
fp114 = fp115 * fp1;
target[((hsi*150+144)*1+lsi)*1] = fp114;
LIBINT2_REALTYPE fp117;
fp117 = fp44 * fp38;
LIBINT2_REALTYPE fp116;
fp116 = fp117 * fp2;
target[((hsi*150+143)*1+lsi)*1] = fp116;
LIBINT2_REALTYPE fp119;
fp119 = fp44 * fp37;
LIBINT2_REALTYPE fp118;
fp118 = fp119 * fp3;
target[((hsi*150+142)*1+lsi)*1] = fp118;
LIBINT2_REALTYPE fp141;
fp141 = fp41 * fp31;
LIBINT2_REALTYPE fp140;
fp140 = fp141 * fp8;
target[((hsi*150+131)*1+lsi)*1] = fp140;
LIBINT2_REALTYPE fp123;
fp123 = fp47 * fp39;
LIBINT2_REALTYPE fp122;
fp122 = fp123 * fp2;
target[((hsi*150+140)*1+lsi)*1] = fp122;
LIBINT2_REALTYPE fp125;
fp125 = fp47 * fp38;
LIBINT2_REALTYPE fp124;
fp124 = fp125 * fp3;
target[((hsi*150+139)*1+lsi)*1] = fp124;
LIBINT2_REALTYPE fp127;
fp127 = fp47 * fp37;
LIBINT2_REALTYPE fp126;
fp126 = fp127 * fp4;
target[((hsi*150+138)*1+lsi)*1] = fp126;
LIBINT2_REALTYPE fp129;
fp129 = fp50 * fp39;
LIBINT2_REALTYPE fp128;
fp128 = fp129 * fp3;
target[((hsi*150+137)*1+lsi)*1] = fp128;
LIBINT2_REALTYPE fp131;
fp131 = fp50 * fp38;
LIBINT2_REALTYPE fp130;
fp130 = fp131 * fp4;
target[((hsi*150+136)*1+lsi)*1] = fp130;
LIBINT2_REALTYPE fp133;
fp133 = fp54 * fp39;
LIBINT2_REALTYPE fp132;
fp132 = fp133 * fp4;
target[((hsi*150+135)*1+lsi)*1] = fp132;
LIBINT2_REALTYPE fp5;
fp5 = 0 + src2[((hsi*24+16)*1+lsi)*1];
LIBINT2_REALTYPE fp135;
fp135 = fp41 * fp34;
LIBINT2_REALTYPE fp134;
fp134 = fp135 * fp5;
target[((hsi*150+134)*1+lsi)*1] = fp134;
LIBINT2_REALTYPE fp137;
fp137 = fp41 * fp33;
LIBINT2_REALTYPE fp136;
fp136 = fp137 * fp6;
target[((hsi*150+133)*1+lsi)*1] = fp136;
LIBINT2_REALTYPE fp139;
fp139 = fp41 * fp32;
LIBINT2_REALTYPE fp138;
fp138 = fp139 * fp7;
target[((hsi*150+132)*1+lsi)*1] = fp138;
LIBINT2_REALTYPE fp121;
fp121 = fp44 * fp36;
LIBINT2_REALTYPE fp120;
fp120 = fp121 * fp4;
target[((hsi*150+141)*1+lsi)*1] = fp120;
LIBINT2_REALTYPE fp22;
fp22 = 0 + src1[((hsi*24+20)*1+lsi)*1];
LIBINT2_REALTYPE fp217;
fp217 = fp47 * fp22;
LIBINT2_REALTYPE fp216;
fp216 = fp217 * fp19;
target[((hsi*150+93)*1+lsi)*1] = fp216;
LIBINT2_REALTYPE fp24;
fp24 = 0 + src1[((hsi*24+18)*1+lsi)*1];
LIBINT2_REALTYPE fp219;
fp219 = fp50 * fp24;
LIBINT2_REALTYPE fp218;
fp218 = fp219 * fp18;
target[((hsi*150+92)*1+lsi)*1] = fp218;
LIBINT2_REALTYPE fp23;
fp23 = 0 + src1[((hsi*24+19)*1+lsi)*1];
LIBINT2_REALTYPE fp221;
fp221 = fp50 * fp23;
LIBINT2_REALTYPE fp220;
fp220 = fp221 * fp19;
target[((hsi*150+91)*1+lsi)*1] = fp220;
LIBINT2_REALTYPE fp223;
fp223 = fp54 * fp24;
LIBINT2_REALTYPE fp222;
fp222 = fp223 * fp19;
target[((hsi*150+90)*1+lsi)*1] = fp222;
LIBINT2_REALTYPE fp225;
fp225 = fp57 * fp39;
LIBINT2_REALTYPE fp224;
fp224 = fp225 * fp5;
target[((hsi*150+89)*1+lsi)*1] = fp224;
LIBINT2_REALTYPE fp227;
fp227 = fp57 * fp38;
LIBINT2_REALTYPE fp226;
fp226 = fp227 * fp6;
target[((hsi*150+88)*1+lsi)*1] = fp226;
LIBINT2_REALTYPE fp229;
fp229 = fp57 * fp37;
LIBINT2_REALTYPE fp228;
fp228 = fp229 * fp7;
target[((hsi*150+87)*1+lsi)*1] = fp228;
LIBINT2_REALTYPE fp231;
fp231 = fp57 * fp36;
LIBINT2_REALTYPE fp230;
fp230 = fp231 * fp8;
target[((hsi*150+86)*1+lsi)*1] = fp230;
LIBINT2_REALTYPE fp177;
fp177 = fp44 * fp28;
LIBINT2_REALTYPE fp176;
fp176 = fp177 * fp12;
target[((hsi*150+113)*1+lsi)*1] = fp176;
LIBINT2_REALTYPE fp235;
fp235 = fp60 * fp39;
LIBINT2_REALTYPE fp234;
fp234 = fp235 * fp6;
target[((hsi*150+84)*1+lsi)*1] = fp234;
LIBINT2_REALTYPE fp237;
fp237 = fp60 * fp38;
LIBINT2_REALTYPE fp236;
fp236 = fp237 * fp7;
target[((hsi*150+83)*1+lsi)*1] = fp236;
LIBINT2_REALTYPE fp239;
fp239 = fp60 * fp37;
LIBINT2_REALTYPE fp238;
fp238 = fp239 * fp8;
target[((hsi*150+82)*1+lsi)*1] = fp238;
LIBINT2_REALTYPE fp241;
fp241 = fp60 * fp36;
LIBINT2_REALTYPE fp240;
fp240 = fp241 * fp9;
target[((hsi*150+81)*1+lsi)*1] = fp240;
LIBINT2_REALTYPE fp243;
fp243 = fp63 * fp39;
LIBINT2_REALTYPE fp242;
fp242 = fp243 * fp7;
target[((hsi*150+80)*1+lsi)*1] = fp242;
LIBINT2_REALTYPE fp245;
fp245 = fp63 * fp38;
LIBINT2_REALTYPE fp244;
fp244 = fp245 * fp8;
target[((hsi*150+79)*1+lsi)*1] = fp244;
LIBINT2_REALTYPE fp247;
fp247 = fp63 * fp37;
LIBINT2_REALTYPE fp246;
fp246 = fp247 * fp9;
target[((hsi*150+78)*1+lsi)*1] = fp246;
LIBINT2_REALTYPE fp249;
fp249 = fp66 * fp39;
LIBINT2_REALTYPE fp248;
fp248 = fp249 * fp8;
target[((hsi*150+77)*1+lsi)*1] = fp248;
LIBINT2_REALTYPE fp251;
fp251 = fp66 * fp38;
LIBINT2_REALTYPE fp250;
fp250 = fp251 * fp9;
target[((hsi*150+76)*1+lsi)*1] = fp250;
LIBINT2_REALTYPE fp233;
fp233 = fp57 * fp35;
LIBINT2_REALTYPE fp232;
fp232 = fp233 * fp9;
target[((hsi*150+85)*1+lsi)*1] = fp232;
LIBINT2_REALTYPE fp179;
fp179 = fp44 * fp27;
LIBINT2_REALTYPE fp178;
fp178 = fp179 * fp13;
target[((hsi*150+112)*1+lsi)*1] = fp178;
LIBINT2_REALTYPE fp181;
fp181 = fp44 * fp26;
LIBINT2_REALTYPE fp180;
fp180 = fp181 * fp14;
target[((hsi*150+111)*1+lsi)*1] = fp180;
LIBINT2_REALTYPE fp183;
fp183 = fp47 * fp29;
LIBINT2_REALTYPE fp182;
fp182 = fp183 * fp12;
target[((hsi*150+110)*1+lsi)*1] = fp182;
LIBINT2_REALTYPE fp185;
fp185 = fp47 * fp28;
LIBINT2_REALTYPE fp184;
fp184 = fp185 * fp13;
target[((hsi*150+109)*1+lsi)*1] = fp184;
LIBINT2_REALTYPE fp187;
fp187 = fp47 * fp27;
LIBINT2_REALTYPE fp186;
fp186 = fp187 * fp14;
target[((hsi*150+108)*1+lsi)*1] = fp186;
LIBINT2_REALTYPE fp189;
fp189 = fp50 * fp29;
LIBINT2_REALTYPE fp188;
fp188 = fp189 * fp13;
target[((hsi*150+107)*1+lsi)*1] = fp188;
LIBINT2_REALTYPE fp191;
fp191 = fp50 * fp28;
LIBINT2_REALTYPE fp190;
fp190 = fp191 * fp14;
target[((hsi*150+106)*1+lsi)*1] = fp190;
LIBINT2_REALTYPE fp193;
fp193 = fp54 * fp29;
LIBINT2_REALTYPE fp192;
fp192 = fp193 * fp14;
target[((hsi*150+105)*1+lsi)*1] = fp192;
LIBINT2_REALTYPE fp215;
fp215 = fp47 * fp23;
LIBINT2_REALTYPE fp214;
fp214 = fp215 * fp18;
target[((hsi*150+94)*1+lsi)*1] = fp214;
LIBINT2_REALTYPE fp197;
fp197 = fp41 * fp23;
LIBINT2_REALTYPE fp196;
fp196 = fp197 * fp16;
target[((hsi*150+103)*1+lsi)*1] = fp196;
LIBINT2_REALTYPE fp199;
fp199 = fp41 * fp22;
LIBINT2_REALTYPE fp198;
fp198 = fp199 * fp17;
target[((hsi*150+102)*1+lsi)*1] = fp198;
LIBINT2_REALTYPE fp21;
fp21 = 0 + src1[((hsi*24+21)*1+lsi)*1];
LIBINT2_REALTYPE fp201;
fp201 = fp41 * fp21;
LIBINT2_REALTYPE fp200;
fp200 = fp201 * fp18;
target[((hsi*150+101)*1+lsi)*1] = fp200;
LIBINT2_REALTYPE fp20;
fp20 = 0 + src1[((hsi*24+22)*1+lsi)*1];
LIBINT2_REALTYPE fp203;
fp203 = fp41 * fp20;
LIBINT2_REALTYPE fp202;
fp202 = fp203 * fp19;
target[((hsi*150+100)*1+lsi)*1] = fp202;
LIBINT2_REALTYPE fp205;
fp205 = fp44 * fp24;
LIBINT2_REALTYPE fp204;
fp204 = fp205 * fp16;
target[((hsi*150+99)*1+lsi)*1] = fp204;
LIBINT2_REALTYPE fp207;
fp207 = fp44 * fp23;
LIBINT2_REALTYPE fp206;
fp206 = fp207 * fp17;
target[((hsi*150+98)*1+lsi)*1] = fp206;
LIBINT2_REALTYPE fp209;
fp209 = fp44 * fp22;
LIBINT2_REALTYPE fp208;
fp208 = fp209 * fp18;
target[((hsi*150+97)*1+lsi)*1] = fp208;
LIBINT2_REALTYPE fp211;
fp211 = fp44 * fp21;
LIBINT2_REALTYPE fp210;
fp210 = fp211 * fp19;
target[((hsi*150+96)*1+lsi)*1] = fp210;
LIBINT2_REALTYPE fp213;
fp213 = fp47 * fp24;
LIBINT2_REALTYPE fp212;
fp212 = fp213 * fp17;
target[((hsi*150+95)*1+lsi)*1] = fp212;
LIBINT2_REALTYPE fp195;
fp195 = fp41 * fp24;
LIBINT2_REALTYPE fp194;
fp194 = fp195 * fp15;
target[((hsi*150+104)*1+lsi)*1] = fp194;
}
}
}
const int hsi = 0;
const int lsi = 0;
const int vi = 0;
/** Number of flops = 404 */
}

#ifdef __cplusplus
};
#endif
