#include "piv_ge_solver.h"
#include <stdlib.h>
#include <gsl/gsl_matrix_double.h> //dolaczamy biblioteke GSL

int piv_ge_solver (matrix_t * eqs){

	if (eqs != NULL) {
		//pivot_ge_in_situ_matrix (eqs);
		gsl_matrix * m = gsl_matrix_alloc (eqs->rn, eqs->cn);
		m->size1 = eqs->rn;
		m->size2 = eqs->cn;
		m->data = eqs->e;

		int i, k;
		int j = 0;
		for (k = 0; k < m->size1 - 1; k++) { /* eliminujemy (zerujemy) kolumnę nr k */
			int piv = k; /* wybór elementu dominującego - maks. z k-tej kol., poniżej diag */
			for (i = k + 1; i < m->size1; i++)
			if (fabs (*(m->data + i * m->size2 + k)) > fabs (*(m->data + piv * m->size2 + k)))
			piv = i;
			// void gsl_matrix_max_index (m,i, j);
			if (piv != k) /* jeśli diag. nta + piv * m->size2 + k)tae jest pivtem - wymień wiersze */
			gsl_matrix_swap_rows (m, piv, k);

			for (i = k + 1; i < m->size1; i++) { /* pętla po kolejnych wierszach poniżej diagonalii k,k */
				double d = *((m->data + i * m->size2 + k)) / *((m->data + piv * m->size2 + k));
				for (j = k; j < m->size2; j++)
				*(m->data + i * m->size2 + j) -= d * *(m->data + k * m->size2 + j);
			}
		}
		eqs->rn = m->size1;
		eqs->cn = m->size2;
		eqs->e = m->data;

		gsl_matrix_free(m);
		if (bs_matrix (eqs) == 0)
		return 0;

		else {
		return 1;
		}
	}
else
return 1;
}
