/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#define MAC_Apply_G_1x2_ops( gamma, sigma, beta, epsilon ) \
{ \
	*(beta)    = *(epsilon) * *(sigma); \
	*(epsilon) = *(epsilon) * *(gamma); \
}

#define MAC_Apply_G_1x2_opd( gamma, sigma, beta, epsilon ) \
{ \
	*(beta)    = *(epsilon) * *(sigma); \
	*(epsilon) = *(epsilon) * *(gamma); \
}

