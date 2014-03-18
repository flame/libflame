
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

