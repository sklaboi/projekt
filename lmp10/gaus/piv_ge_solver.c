#include "piv_ge_solver.h"
#include <stdlib.h>
#include <gsl/gsl_matrix_double.h> //dolaczamy biblioteke GSL

int piv_ge_solver (matrix_t * eqs){

	if (eqs != NULL) {
		gsl_matrix * m = gsl_matrix_alloc (eqs->rn, eqs->cn); //Ta funkcja gsl_matrix_alloc(size_t n1, size_t n2) tworze macierz w rozmierze n1xn2, zwracajac wskaznik do nowo zainicjowanej struktury macierzy
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
			gsl_matrix_swap_rows (m, piv, k);//P.W.: int gsl_matrix_swap_rows(gsl_matrix * M, size_t i, size_t j):Funkcja ta wymienia i-tego i j-tej wiersze matrycy w mejscu M;

			for (i = k + 1; i < m->size1; i++) { /* pętla po kolejnych wierszach poniżej diagonalii k,k */
				double d = *((m->data + i * m->size2 + k)) / *((m->data + piv * m->size2 + k));
				for (j = k; j < m->size2; j++)
				*(m->data + i * m->size2 + j) -= d * *(m->data + k * m->size2 + j);
			}
		}
		eqs->rn = m->size1;
		eqs->cn = m->size2;
		eqs->e = m->data;

		gsl_matrix_free(m);//Ta funkcja zwalnia wczesniej przydzielona macierz m; U nas matryca zostala utworzona za pomoca gsl_matrix_alloc, to nastepnie "blok" bazowej matrycy beda tez zwalniane
		if (bs_matrix (eqs) == 0)
		return 0;

		else {
		return 1;
		}
	}
else
return 1;
}
