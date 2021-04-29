#include <stdlib.h>
#include <cmath>

/*************************************************************************
 ROCKFLOW - Funktion: Gauss

 Aufgabe:
   Gleichungsloeser:
   LR-Faktorisierung von matrix nach dem Gausschen Algorithmus mit
   partieller Pivotisierung; Berechnung der Loesung aus LR-Faktorisierung
   Der Ergebnisvektor wird erst auf die rechte Seite gespeichert, und
   hinterher umkopiert. Hier kann noch Speicherplatz gespart werden.
   Die alten Inhalte von matrix und vecb werden zerstoert!
   Quelle: Schwetlick/Kretzschmar, Numerische Verfahren 1991

 Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E double *matrix: Linke Seite Gleichungssystem
   E double *vecb: Rechte Seite Gleichungssystem
   X double *vecx: Ergebnisvektor, Speicher muss bereits reserviert sein
   E int g: Dimension des Gleichungssystems

 Ergebnis:
   - void -

 Programmaenderungen:
   12/1994     MSR        Erste Version

*************************************************************************/
void Gauss(double *matrix, double *vecb, double *vecx, int g)
{
    /* Matrizen sind in C zeilenweise abgespeichert:
       matrix[i][j] -> matrix[i*g+j] */
    static int *s;
    static double z, hilf;
    register int k, i, j, sk;
    //s = (int *) Malloc(sizeof(int) * (g - 1));
    s = new int[g-1];
    /* LR-Faktorisierung */
    for (k = 0; k < (g - 1); k++) {
        /* Pivotsuche */
        z = 0.0;
        sk = 0;
        for (i = k; i < g; i++) {
            hilf = fabs(matrix[i * g + k]);     /* matrix[i][k] */
            if (hilf > z) {
                z = hilf;
                sk = i;
            }
        }
        s[k] = sk;
        /* evtl. Zeilen vertauschen */
        if (sk > k) {
            for (j = 0; j < g; j++) {
                z = matrix[k * g + j];  /* matrix[k][j] */
                matrix[k * g + j] = matrix[sk * g + j];         /* matrix[k][j], matrix[sk][j] */
                matrix[sk * g + j] = z;         /* matrix[sk][j] */
            }
        }
        /* Berechnung der Eliminationskoeffizienten */
        for (i = (k + 1); i < g; i++) {
            matrix[i * g + k] /= matrix[k * g + k];     /* matrix[i][k], matrix[k][k] */
        }
        /* Spaltenweise Berechnung der neuen Restmatrix */
        for (j = (k + 1); j < g; j++) {
            for (i = (k + 1); i < g; i++) {
                matrix[i * g + j] -= (matrix[i * g + k] * matrix[k * g + j]);
                /* matrix[i][j], matrix[i][k], matrix[k][j] */
            }
        }
    }
    /* Loesung berechnen */
    /* vecb transformieren */
    for (k = 0; k < (g - 1); k++) {
        sk = s[k];
        if (sk > k) {
            z = vecb[k];
            vecb[k] = vecb[sk];
            vecb[sk] = z;
        }
    }
    for (k = 1; k < g; k++) {
        for (j = 0; j < k; j++) {
            vecb[k] -= (matrix[k * g + j] * vecb[j]);   /* matrix[k][j] */
        }
    }
    /* vecx berechnen */
    for (k = (g - 1); k >= 0; k--) {
        for (j = (k + 1); j < g; j++) {
            vecb[k] -= (matrix[k * g + j] * vecb[j]);   /* matrix[k][j] */
        }
        vecb[k] /= matrix[k * g + k];   /* matrix[k][k] */
    }
    /* Umspeichern des Ergebnisses von vecb nach vecx */
    for (k = 0; k < g; k++)
        vecx[k] = vecb[k];
    /* Speicher freigeben */
    //s = (int *) Free(s);
    delete [] s;
}


