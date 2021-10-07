/*****************************************************************************
 * This file is part of Kvazaar HEVC encoder.
 *
 * Copyright (c) 2021, Tampere University, ITU/ISO/IEC, project contributors
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 * 
 * * Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 * 
 * * Redistributions in binary form must reproduce the above copyright notice, this
 *   list of conditions and the following disclaimer in the documentation and/or
 *   other materials provided with the distribution.
 * 
 * * Neither the name of the Tampere University or ITU/ISO/IEC nor the names of its
 *   contributors may be used to endorse or promote products derived from
 *   this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * INCLUDING NEGLIGENCE OR OTHERWISE ARISING IN ANY WAY OUT OF THE USE OF THIS
 ****************************************************************************/

#include "ml_classifier_intra_depth_pred.h"


int tree_predict_merge_depth_1(features_s* p_features, double* p_nb_iter, double* p_nb_bad)
{
	if (p_features->merge_variance <= 140.3129)
	{
		if (p_features->var_of_sub_var <= 569.6553)
		{
			if (p_features->merge_variance <= 20.8854)
			{
				*p_nb_iter = 19428.0;
				*p_nb_bad = 1740.0;
				return -1.0000;
			}
			else if (p_features->sub_variance_0 <= 9.1015)
			{
				if (p_features->merge_variance <= 39.132)
				{
					*p_nb_iter = 1166.0;
					*p_nb_bad = 358.0;
					return -1.0000;
				}
				else {
					*p_nb_iter = 1049.0;
					*p_nb_bad = 392.0;
					return 1.0000;
				}
			}
			else {
				*p_nb_iter = 9371.0;
				*p_nb_bad = 1805.0;
				return -1.0000;
			}
		}
		else if (p_features->sub_variance_2 <= 23.3193)
		{
			*p_nb_iter = 1059.0;
			*p_nb_bad = 329.0;
			return 1.0000;
		}
		else if (p_features->sub_variance_1 <= 30.7348)
		{
			*p_nb_iter = 1042.0;
			*p_nb_bad = 395.0;
			return 1.0000;
		}
		else {
			*p_nb_iter = 1756.0;
			*p_nb_bad = 588.0;
			return -1.0000;
		}
	}
	else if (p_features->merge_variance <= 857.8047)
	{
		if (p_features->var_of_sub_var <= 66593.5553)
		{
			if (p_features->sub_variance_0 <= 12.1697)
			{
				*p_nb_iter = 2006.0;
				*p_nb_bad = 374.0;
				return 1.0000;
			}
			else if (p_features->neigh_variance_C <= 646.8204)
			{
				if (p_features->neigh_variance_A <= 664.7609)
				{
					if (p_features->neigh_variance_B <= 571.2004)
					{
						if (p_features->var_of_sub_mean <= 4.1069)
						{
							*p_nb_iter = 1208.0;
							*p_nb_bad = 399.0;
							return 1.0000;
						}
						else if (p_features->var_of_sub_var <= 11832.6635)
						{
							*p_nb_iter = 8701.0;
							*p_nb_bad = 3037.0;
							return -1.0000;
						}
						else if (p_features->neigh_variance_A <= 142.298)
						{
							*p_nb_iter = 1025.0;
							*p_nb_bad = 290.0;
							return 1.0000;
						}
						else if (p_features->variance <= 394.4839)
						{
							*p_nb_iter = 1156.0;
							*p_nb_bad = 489.0;
							return 1.0000;
						}
						else {
							*p_nb_iter = 1150.0;
							*p_nb_bad = 503.0;
							return -1.0000;
						}
					}
					else {
						*p_nb_iter = 1777.0;
						*p_nb_bad = 558.0;
						return 1.0000;
					}
				}
				else {
					*p_nb_iter = 1587.0;
					*p_nb_bad = 411.0;
					return 1.0000;
				}
			}
			else {
				*p_nb_iter = 1980.0;
				*p_nb_bad = 474.0;
				return 1.0000;
			}
		}
		else {
			*p_nb_iter = 3613.0;
			*p_nb_bad = 475.0;
			return 1.0000;
		}
	}
	else {
		*p_nb_iter = 20926.0;
		*p_nb_bad = 1873.0;
		return 1.0000;
	}
}



int tree_predict_merge_depth_2(features_s* p_features, double* p_nb_iter, double* p_nb_bad)
{
	if (p_features->merge_variance <= 119.4611)
	{
		if (p_features->var_of_sub_var <= 1078.0638)
		{
			if (p_features->neigh_variance_B <= 70.2189)
			{
				*p_nb_iter = 29253.0;
				*p_nb_bad = 3837.0;
				return -1.0000;
			}
			else if (p_features->variance <= 20.8711)
			{
				*p_nb_iter = 1292.0;
				*p_nb_bad = 458.0;
				return 2.0000;
			}
			else {
				*p_nb_iter = 1707.0;
				*p_nb_bad = 399.0;
				return -1.0000;
			}
		}
		else if (p_features->var_of_sub_var <= 3300.4034)
		{
			*p_nb_iter = 1554.0;
			*p_nb_bad = 675.0;
			return -1.0000;
		}
		else {
			*p_nb_iter = 1540.0;
			*p_nb_bad = 429.0;
			return 2.0000;
		}
	}
	else if (p_features->merge_variance <= 696.1989)
	{
		if (p_features->var_of_sub_var <= 31803.3242)
		{
			if (p_features->sub_variance_2 <= 10.3845)
			{
				*p_nb_iter = 3473.0;
				*p_nb_bad = 768.0;
				return 2.0000;
			}
			else if (p_features->neigh_variance_C <= 571.5329)
			{
				if (p_features->neigh_variance_B <= 492.8159)
				{
					if (p_features->neigh_variance_B <= 38.9672)
					{
						*p_nb_iter = 1887.0;
						*p_nb_bad = 588.0;
						return 2.0000;
					}
					else if (p_features->neigh_variance_A <= 380.5927)
					{
						if (p_features->sub_variance_1 <= 19.9678)
						{
							*p_nb_iter = 1686.0;
							*p_nb_bad = 721.0;
							return 2.0000;
						}
						else if (p_features->neigh_variance_A <= 66.6749)
						{
							*p_nb_iter = 1440.0;
							*p_nb_bad = 631.0;
							return 2.0000;
						}
						else {
							*p_nb_iter = 5772.0;
							*p_nb_bad = 2031.0;
							return -1.0000;
						}
					}
					else {
						*p_nb_iter = 1791.0;
						*p_nb_bad = 619.0;
						return 2.0000;
					}
				}
				else {
					*p_nb_iter = 1624.0;
					*p_nb_bad = 494.0;
					return 2.0000;
				}
			}
			else {
				*p_nb_iter = 1298.0;
				*p_nb_bad = 312.0;
				return 2.0000;
			}
		}
		else {
			*p_nb_iter = 4577.0;
			*p_nb_bad = 892.0;
			return 2.0000;
		}
	}
	else {
		*p_nb_iter = 21106.0;
		*p_nb_bad = 2744.0;
		return 2.0000;
	}
}



int tree_predict_merge_depth_3(features_s* p_features, double* p_nb_iter, double* p_nb_bad)
{
	if (p_features->merge_variance <= 80.1487)
	{
		if (p_features->neigh_variance_C <= 83.7148)
		{
			*p_nb_iter = 29806.0;
			*p_nb_bad = 3603.0;
			return -1.0000;
		}
		else {
			*p_nb_iter = 1003.0;
			*p_nb_bad = 421.0;
			return 3.0000;
		}
	}
	else if (p_features->merge_variance <= 351.8138)
	{
		if (p_features->neigh_variance_C <= 255.4236)
		{
			if (p_features->neigh_variance_B <= 260.5349)
			{
				if (p_features->var_of_sub_var <= 6381.513)
				{
					if (p_features->neigh_variance_A <= 244.2556)
					{
						if (p_features->sub_variance_0 <= 4.75)
						{
							*p_nb_iter = 1290.0;
							*p_nb_bad = 525.0;
							return 3.0000;
						}
						else if (p_features->neigh_variance_B <= 16.9287)
						{
							*p_nb_iter = 1045.0;
							*p_nb_bad = 499.0;
							return 3.0000;
						}
						else {
							*p_nb_iter = 6901.0;
							*p_nb_bad = 2494.0;
							return -1.0000;
						}
					}
					else {
						*p_nb_iter = 1332.0;
						*p_nb_bad = 408.0;
						return 3.0000;
					}
				}
				else {
					*p_nb_iter = 2929.0;
					*p_nb_bad = 842.0;
					return 3.0000;
				}
			}
			else {
				*p_nb_iter = 2239.0;
				*p_nb_bad = 572.0;
				return 3.0000;
			}
		}
		else {
			*p_nb_iter = 2777.0;
			*p_nb_bad = 714.0;
			return 3.0000;
		}
	}
	else {
		*p_nb_iter = 30678.0;
		*p_nb_bad = 5409.0;
		return 3.0000;
	}
}



int tree_predict_merge_depth_4(features_s* p_features, double* p_nb_iter, double* p_nb_bad)
{
	if (p_features->neigh_variance_C <= 240.2773)
	{
		if (p_features->neigh_variance_B <= 227.5898)
		{
			if (p_features->neigh_variance_A <= 195.4844)
			{
				if (p_features->variance <= 203.3086)
				{
					if (p_features->qp <= 32)
					{
						if (p_features->neigh_variance_C <= 102.2344)
						{
							if (p_features->neigh_variance_B <= 116.4961)
							{
								if (p_features->variance <= 89.4023)
								{
									*p_nb_iter = 27398.0;
									*p_nb_bad = 4665.0;
									return -1.0000;
								}
								else {
									*p_nb_iter = 1676.0;
									*p_nb_bad = 795.0;
									return 4.0000;
								}
							}
							else {
								*p_nb_iter = 1405.0;
								*p_nb_bad = 566.0;
								return 4.0000;
							}
						}
						else {
							*p_nb_iter = 2827.0;
							*p_nb_bad = 1173.0;
							return 4.0000;
						}
					}
					else {
						*p_nb_iter = 8871.0;
						*p_nb_bad = 822.0;
						return -1.0000;
					}
				}
				else {
					*p_nb_iter = 3162.0;
					*p_nb_bad = 718.0;
					return 4.0000;
				}
			}
			else {
				*p_nb_iter = 6154.0;
				*p_nb_bad = 1397.0;
				return 4.0000;
			}
		}
		else {
			*p_nb_iter = 9385.0;
			*p_nb_bad = 1609.0;
			return 4.0000;
		}
	}
	else {
		*p_nb_iter = 19122.0;
		*p_nb_bad = 2960.0;
		return 4.0000;
	}
}



int tree_predict_split_depth_0(features_s* p_features, double* p_nb_iter, double* p_nb_bad)
{
	if (p_features->var_of_sub_var <= 12754.7856)
	{
		if (p_features->var_of_sub_var <= 137.9034)
		{
			*p_nb_iter = 25155.0;
			*p_nb_bad = 2959.0;
			return 0.0000;
		}
		else if (p_features->sub_variance_2 <= 13.2892)
		{
			*p_nb_iter = 1080.0;
			*p_nb_bad = 383.0;
			return -1.0000;
		}
		else if (p_features->variance <= 564.1738)
		{
			if (p_features->var_of_sub_var <= 1185.4728)
			{
				*p_nb_iter = 6067.0;
				*p_nb_bad = 1699.0;
				return 0.0000;
			}
			else if (p_features->var_of_sub_mean <= 46.2388)
			{
				if (p_features->sub_variance_0 <= 46.8708)
				{
					*p_nb_iter = 1088.0;
					*p_nb_bad = 377.0;
					return -1.0000;
				}
				else if (p_features->sub_variance_3 <= 61.4213)
				{
					*p_nb_iter = 1183.0;
					*p_nb_bad = 498.0;
					return -1.0000;
				}
				else {
					*p_nb_iter = 3416.0;
					*p_nb_bad = 1373.0;
					return 0.0000;
				}
			}
			else {
				*p_nb_iter = 3769.0;
				*p_nb_bad = 1093.0;
				return 0.0000;
			}
		}
		else {
			*p_nb_iter = 1036.0;
			*p_nb_bad = 434.0;
			return -1.0000;
		}
	}
	else if (p_features->var_of_sub_var <= 98333.8279)
	{
		if (p_features->variance <= 987.2333)
		{
			if (p_features->var_of_sub_var <= 37261.2896)
			{
				if (p_features->variance <= 238.2248)
				{
					*p_nb_iter = 1323.0;
					*p_nb_bad = 301.0;
					return -1.0000;
				}
				else if (p_features->var_of_sub_var <= 17347.3971)
				{
					*p_nb_iter = 1215.0;
					*p_nb_bad = 550.0;
					return 0.0000;
				}
				else if (p_features->qp <= 22)
				{
					*p_nb_iter = 1000.0;
					*p_nb_bad = 493.0;
					return 0.0000;
				}
				else {
					*p_nb_iter = 2640.0;
					*p_nb_bad = 1121.0;
					return -1.0000;
				}
			}
			else {
				*p_nb_iter = 5188.0;
				*p_nb_bad = 1248.0;
				return -1.0000;
			}
		}
		else {
			*p_nb_iter = 2323.0;
			*p_nb_bad = 274.0;
			return -1.0000;
		}
	}
	else {
		*p_nb_iter = 21357.0;
		*p_nb_bad = 1829.0;
		return -1.0000;
	}
}


int tree_predict_split_depth_1(features_s* p_features, double* p_nb_iter, double* p_nb_bad)
{
	if (p_features->var_of_sub_var <= 1138.9473)
	{
		*p_nb_iter = 32445.0;
		*p_nb_bad = 4580.0;
		return 1.0000;
	}
	else if (p_features->var_of_sub_var <= 27289.2117)
	{
		if (p_features->sub_variance_1 <= 12.0603)
		{
			*p_nb_iter = 1900.0;
			*p_nb_bad = 401.0;
			return -1.0000;
		}
		else if (p_features->var_of_sub_var <= 5841.4773)
		{
			if (p_features->variance <= 72.4175)
			{
				*p_nb_iter = 1000.0;
				*p_nb_bad = 356.0;
				return -1.0000;
			}
			else if (p_features->neigh_variance_A <= 633.8163)
			{
				*p_nb_iter = 5279.0;
				*p_nb_bad = 1961.0;
				return 1.0000;
			}
			else {
				*p_nb_iter = 1176.0;
				*p_nb_bad = 527.0;
				return -1.0000;
			}
		}
		else if (p_features->sub_variance_0 <= 38.3035)
		{
			*p_nb_iter = 1251.0;
			*p_nb_bad = 293.0;
			return -1.0000;
		}
		else if (p_features->neigh_variance_B <= 664.9494)
		{
			if (p_features->sub_variance_3 <= 45.8181)
			{
				*p_nb_iter = 1276.0;
				*p_nb_bad = 471.0;
				return -1.0000;
			}
			else if (p_features->sub_variance_3 <= 404.3086)
			{
				if (p_features->sub_variance_1 <= 99.8715)
				{
					*p_nb_iter = 1005.0;
					*p_nb_bad = 435.0;
					return -1.0000;
				}
				else if (p_features->sub_variance_0 <= 282.3064)
				{
					*p_nb_iter = 1370.0;
					*p_nb_bad = 539.0;
					return 1.0000;
				}
				else {
					*p_nb_iter = 1013.0;
					*p_nb_bad = 495.0;
					return -1.0000;
				}
			}
			else {
				*p_nb_iter = 1000.0;
				*p_nb_bad = 379.0;
				return -1.0000;
			}
		}
		else {
			*p_nb_iter = 2270.0;
			*p_nb_bad = 679.0;
			return -1.0000;
		}
	}
	else {
		*p_nb_iter = 29015.0;
		*p_nb_bad = 3950.0;
		return -1.0000;
	}
}


int tree_predict_split_depth_2(features_s* p_features, double* p_nb_iter, double* p_nb_bad)
{
	if (p_features->var_of_sub_var <= 2597.4529)
	{
		if (p_features->var_of_sub_var <= 146.7734)
		{
			*p_nb_iter = 23216.0;
			*p_nb_bad = 1560.0;
			return 2.0000;
		}
		else if (p_features->merge_variance <= 259.6952)
		{
			*p_nb_iter = 7470.0;
			*p_nb_bad = 1902.0;
			return 2.0000;
		}
		else if (p_features->qp <= 27)
		{
			if (p_features->variance <= 73.9929)
			{
				*p_nb_iter = 1138.0;
				*p_nb_bad = 486.0;
				return -1.0000;
			}
			else {
				*p_nb_iter = 1619.0;
				*p_nb_bad = 716.0;
				return 2.0000;
			}
		}
		else {
			*p_nb_iter = 2425.0;
			*p_nb_bad = 861.0;
			return 2.0000;
		}
	}
	else if (p_features->var_of_sub_var <= 60850.5208)
	{
		if (p_features->var_of_sub_var <= 10144.602)
		{
			if (p_features->neigh_variance_C <= 926.8972)
			{
				if (p_features->sub_variance_0 <= 26.6006)
				{
					*p_nb_iter = 1796.0;
					*p_nb_bad = 586.0;
					return -1.0000;
				}
				else if (p_features->neigh_variance_A <= 493.5849)
				{
					if (p_features->neigh_variance_A <= 72.9516)
					{
						*p_nb_iter = 1326.0;
						*p_nb_bad = 557.0;
						return -1.0000;
					}
					else if (p_features->variance <= 156.4014)
					{
						*p_nb_iter = 1210.0;
						*p_nb_bad = 563.0;
						return -1.0000;
					}
					else {
						*p_nb_iter = 1920.0;
						*p_nb_bad = 817.0;
						return 2.0000;
					}
				}
				else {
					*p_nb_iter = 1106.0;
					*p_nb_bad = 437.0;
					return -1.0000;
				}
			}
			else {
				*p_nb_iter = 1001.0;
				*p_nb_bad = 278.0;
				return -1.0000;
			}
		}
		else {
			*p_nb_iter = 13068.0;
			*p_nb_bad = 3612.0;
			return -1.0000;
		}
	}
	else {
		*p_nb_iter = 22705.0;
		*p_nb_bad = 2687.0;
		return -1.0000;
	}
}



int tree_predict_split_depth_3(features_s* p_features, double* p_nb_iter, double* p_nb_bad)
{
	if (p_features->var_of_sub_var <= 818.5173)
	{
		if (p_features->merge_variance <= 62.7641)
		{
			*p_nb_iter = 20568.0;
			*p_nb_bad = 767.0;
			return 3.0000;
		}
		else if (p_features->qp <= 27)
		{
			if (p_features->variance <= 9.4219)
			{
				*p_nb_iter = 1255.0;
				*p_nb_bad = 206.0;
				return 3.0000;
			}
			else if (p_features->merge_variance <= 375.2185)
			{
				*p_nb_iter = 3999.0;
				*p_nb_bad = 1321.0;
				return 3.0000;
			}
			else {
				*p_nb_iter = 1786.0;
				*p_nb_bad = 817.0;
				return -1.0000;
			}
		}
		else {
			*p_nb_iter = 5286.0;
			*p_nb_bad = 737.0;
			return 3.0000;
		}
	}
	else if (p_features->var_of_sub_var <= 37332.3018)
	{
		if (p_features->var_of_sub_var <= 7585.0282)
		{
			if (p_features->qp <= 32)
			{
				if (p_features->neigh_variance_C <= 330.2178)
				{
					if (p_features->sub_variance_0 <= 8.5273)
					{
						*p_nb_iter = 1114.0;
						*p_nb_bad = 346.0;
						return -1.0000;
					}
					else if (p_features->neigh_variance_B <= 221.5469)
					{
						if (p_features->var_of_sub_var <= 1989.7928)
						{
							*p_nb_iter = 1539.0;
							*p_nb_bad = 606.0;
							return 3.0000;
						}
						else if (p_features->variance <= 155.5974)
						{
							*p_nb_iter = 1298.0;
							*p_nb_bad = 634.0;
							return 3.0000;
						}
						else {
							*p_nb_iter = 1076.0;
							*p_nb_bad = 456.0;
							return -1.0000;
						}
					}
					else {
						*p_nb_iter = 1644.0;
						*p_nb_bad = 639.0;
						return -1.0000;
					}
				}
				else {
					*p_nb_iter = 2401.0;
					*p_nb_bad = 713.0;
					return -1.0000;
				}
			}
			else if (p_features->merge_variance <= 281.9509)
			{
				*p_nb_iter = 1020.0;
				*p_nb_bad = 262.0;
				return 3.0000;
			}
			else {
				*p_nb_iter = 1278.0;
				*p_nb_bad = 594.0;
				return -1.0000;
			}
		}
		else {
			*p_nb_iter = 10507.0;
			*p_nb_bad = 2943.0;
			return -1.0000;
		}
	}
	else {
		*p_nb_iter = 25229.0;
		*p_nb_bad = 3060.0;
		return -1.0000;
	}
}
