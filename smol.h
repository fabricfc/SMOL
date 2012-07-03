#include <stdio.h>
#include <cstdlib>

/**
 * SMOL
 * Some Matrix Operations Library (creates, multiplicates, transposes, inverse, destroy)
 * Author: Fabricio F. Carvalho
 * License: GNU
 * 
 * thanks for Jason Yu-Tseh (inversion code)
 * http://chi3x10.wordpress.com/2008/05/28/calculate-matrix-inversion-in-c/
 */

#ifndef SMOL_H
#define	SMOL_H

using namespace std;

template <class T>
class SMOL {
private:

    // calculate the cofactor of element (row,col)
    void GetMinor(T **src, T **dest, int row, int col, int order);

public:
    /**
     * matrix multiplication 
     * multiplica a matriz
     */
    T **matrizMult(T **mat1, int lin1, int col1, T **mat2, int lin2, int col2);
    /**
     * matrix multiplication 
     * multiplica a matriz
     * @param matResult Result of the operation
     */
    void matrizMult(T **mat1, int lin1, int col1, T **mat2, int lin2, int col2, T **matResult);
    /**
     * matrix sum 
     * soma de matriz
     */
    T **matrizSum(T **mat1, int lin1, int col1, T **mat2, int lin2, int col2);
    /**
     * matrix sum 
     * soma de matriz
     * @param matResult Result of the operation
     */
    void matrizSum(T **mat1, int lin1, int col1, T **mat2, int lin2, int col2, T **matResult);
    /**
     * matrix subtraction
     * subtração of matriz
     */
    T **matrizSub(T **mat1, int lin1, int col1, T **mat2, int lin2, int col2);
    /**
     * matrix subtraction
     * subtração of matriz
     * @param matResult Result of the operation
     */
    void matrizSub(T **mat1, int lin1, int col1, T **mat2, int lin2, int col2, T **matResult);
    /**
     * matrix transposing
     * transposta a matriz
     */
    T **matrizTrans(T **mat1, int lin1, int col1);
    /**
     * matrix transposing
     * transposta a matriz
     * @param matResult Result of the operation
     */    
    void matrizTrans(T **mat1, int lin1, int col1, T **matResult);
    /**
     * matrix printing
     * imprime a matriz
     */
    void printMatriz(T **mat1, int lines, int cols, char *stringPrintf);
    /**
     * matrix creation
     * cria a matriz
     */
    T **criarMatriz(int lines, int cols);
    /**
     * matrix deallocation 
     * desaloca a matriz
     */
    void desalocaMatriz(T **matriz, int cols);

    // Calcula o determinando recursivamente
    // Calculate the determinant recursively.
    float calculaDeterminante(T **mat, int order);

    // matrix inversioon
    // the result is put in matResult
    void inversaoMatriz(T **A, int order, T **matResult);
};

// calculate the cofactor of element (row,col)

template <class T>
void SMOL<T>::GetMinor (T **src, T **dest, int row, int col, int order) {
    // indicate which col and row is being copied to dest
    int colCount = 0, rowCount = 0;

    for (int i = 0; i < order; i++) {
        if (i != row) {
            colCount = 0;
            for (int j = 0; j < order; j++) {
                // when j is not the element
                if (j != col) {
                    dest[rowCount][colCount] = src[i][j];
                    colCount++;
                }
            }
            rowCount++;
        }
    }
}

// Calculate the determinant recursively.

template <class T>
float SMOL<T>::calculaDeterminante (T **mat, int order) {
    // order must be >= 0
    // stop the recursion when matrix is a single element
    if (order == 1)
        return mat[0][0];

    // the determinant value
    T det = 0;

    // allocate the cofactor matrix
    T **minor;
    minor = new T*[order - 1];
    for (int i = 0; i < order - 1; i++)
        minor[i] = new T[order - 1];

    for (int i = 0; i < order; i++) {
        // get minor of element (0,i)
        GetMinor(mat, minor, 0, i, order);
        // the recusion is here!

        det += (i % 2 == 1 ? -1.0 : 1.0) * mat[0][i] * calculaDeterminante(minor, order - 1);
        //det += pow( -1.0, i ) * mat[0][i] * calculaDeterminante( minor,order-1 );
    }

    // release memory
    for (int i = 0; i < order - 1; i++)
        delete [] minor[i];
    delete [] minor;

    return det;
}


// matrix inversioon
// the result is put in matResult

template <class T>
void SMOL<T>::inversaoMatriz (T **A, int order, T **matResult) {
    // get the determinant of a
    T det = 1.0 / calculaDeterminante(A, order);

    // memory allocation
    T *temp = new T[(order - 1)*(order - 1)];
    T **minor = new T*[order - 1];
    for (int i = 0; i < order - 1; i++)
        minor[i] = temp + (i * (order - 1));

    for (int j = 0; j < order; j++) {
        for (int i = 0; i < order; i++) {
            // get the co-factor (matrix) of A(j,i)
            GetMinor(A, minor, j, i, order);
            matResult[i][j] = det * calculaDeterminante(minor, order - 1);
            if ((i + j) % 2 == 1)
                matResult[i][j] = -matResult[i][j];
        }
    }

    // release memory
    //delete [] minor[0];
    delete [] temp;
    delete [] minor;
}

template <class T>
void SMOL<T>::matrizMult (T **mat1, int lin1, int col1, T **mat2, int lin2, int col2, T **matResult) {
    if (col1 != lin2)
        //multiplicacao incompatível
        return;
    int i, j, k;
    T sum;
    for (i = 0; i < lin1; i++)
        for (j = 0; j < col2; j++) {
            sum = 0;
            for (k = 0; k < col1; k++)
                sum = sum + mat1[i][k] * mat2[k][j];
            matResult[i][j] = sum;
        }
}

template <class T>
T **SMOL<T>::matrizMult (T **mat1, int lin1, int col1, T **mat2, int lin2, int col2) {
    if (col1 != lin2)
        //multiplicacao incompatível
        return NULL;
    T **resultado = criarMatriz(lin1, col2);
    matrizMult(mat1, lin1, col1, mat2, lin2, col2, resultado);
    return resultado;
}

template <class T>
void SMOL<T>::matrizSum (T **mat1, int lin1, int col1, T **mat2, int lin2, int col2, T **matResult) {
    if (col1 != col2 || lin1 != lin2)
        //soma incompatível
        return;
    int i, j, k;
    for (i = 0; i < lin1; i++)
        for (j = 0; j < col1; j++) {
            matResult[i][j] = mat1[i][j] + mat2[i][j];
        }
}

template <class T>
T **SMOL<T>::matrizSum (T **mat1, int lin1, int col1, T **mat2, int lin2, int col2) {
    if (col1 != col2 || lin1 != lin2)
        //soma incompatível
        return NULL;
    T **resultado = criarMatriz(lin1, col1);
    matrizSum(mat1, lin1, col1, mat2, lin2, col2, resultado);
    return resultado;
}

template <class T>
void SMOL<T>::matrizSub (T **mat1, int lin1, int col1, T **mat2, int lin2, int col2, T **matResult)
{
    if (col1 != col2 || lin1 != lin2)
        //soma incompatível
        return;
    int i, j, k;
    for (i = 0; i < lin1; i++)
        for (j = 0; j < col1; j++) {
            matResult[i][j] = mat1[i][j] - mat2[i][j];
        }
}

template <class T>
T **SMOL<T>::matrizSub (T **mat1, int lin1, int col1, T **mat2, int lin2, int col2) {
    if (col1 != col2 || lin1 != lin2)
        //soma incompatível
        return NULL;
    T **resultado = criarMatriz(lin1, col1);
    matrizSub(mat1, lin1, col1, mat2, lin2, col2, resultado);
    return resultado;
}

template <class T>
void SMOL<T>::matrizTrans (T **mat1, int lin1, int col1, T **matResult) {
    int i, j;
    for (i = 0; i < lin1; i++)
        for (j = 0; j < col1; j++) {
            matResult[j][i] = mat1[i][j];
        }
}

template <class T>
T **SMOL<T>::matrizTrans (T **mat1, int lin1, int col1) {
    T **resultado = criarMatriz(col1, lin1);
    matrizTrans(mat1, lin1, col1, resultado);
    return resultado;
}

template <class T>
void SMOL<T>::printMatriz(T **matriz, int lines, int cols, char *stringPrintf) {
    for (int i = 0; i < lines; i++) {
        for (int j = 0; j < cols; j++) {
            printf(stringPrintf, matriz[i][j]);
        }
        printf("\n");
    }
}

template <class T>
T **SMOL<T>::criarMatriz(int lines, int cols) {
    T **matriz;
    matriz = (T**) malloc(lines * sizeof (T*));
    for (int i = 0; i < lines; i++) {
        matriz[i] = (T*) malloc(cols * sizeof (T));
    }
    return matriz;
}

template <class T>
void SMOL<T>::desalocaMatriz(T **matriz, int cols) {
    //desalocando a matriz
    for (int i = 0; i < cols; i++) {
        free(matriz[i]);
    }
    free(matriz);
}

#endif	/* SMOL_H */