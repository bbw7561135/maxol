#ifndef VECTOR_H_
#define VECTOR_H_

#include <math.h>

struct vector {
	double _0;
	double _1;
	double _2;
};

static inline double inner_product(struct vector a, struct vector b)
{
	return a._0 * b._0 + a._1 * b._1 + a._2 * b._2;
}

static inline struct vector exterior_product(struct vector a, struct vector b)
{
	return {a._1 * b._2 - a._2 * b._1,
			a._2 * b._0 - a._0 * b._2,
			a._0 * b._1 - a._1 * b._0};
}

static inline double vector_length(struct vector a)
{
	return sqrt(inner_product(a, a));
}

static inline struct vector vector_add(struct vector a, struct vector b)
{
	return {a._0 + b._0, a._1 + b._1, a._2 + b._2};
}

static inline struct vector vector_multiple(struct vector a, double k)
{
	return {a._0 * k, a._1 * k, a._2 * k};
}

#endif /* VECTOR_H_ */
