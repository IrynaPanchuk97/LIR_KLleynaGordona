using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DirectKleinGordonInteriorNeumann
{
    public static class MatrixHelper
    {
        public static double[,] MatrixIdentity(int n, double index = 1)
        {
            double[,] result = new double[n, n];
            for (int i = 0; i < n; ++i)
                result[i, i] = 1.0 * index;

            return result;
        }
        public static double[] MatrixVector2(double[,] A, double[] B)
        {
            double[] res = new double[B.GetLength(0)];
            for (int row = 0; row < A.GetLength(0); row++)
            {
                for (int col = 0; col < A.GetLength(1); col++)
                {
                    res[col] += A[row, col] * B[row];
                }
            }
            return res;
        }
        public static double[] MatrixVector(double[,] A, double[] B)
        {
            int n = A.GetLength(0);
            int m = A.GetLength(1);
            double[] res = new double[m];
            for (int row = 0; row < n; row++)
            {
                for (int col = 0; col < m; col++)
                {
                    res[col] += A[row, col] * B[col];
                }
            }
            return res;
        }

        public static double[,] Multiplication(double[,] a, double[,] b)
        {
            if (a.GetLength(1) != b.GetLength(0)) throw new Exception("Матриці не можна перемножити");
            double[,] r = new double[a.GetLength(0), b.GetLength(1)];
            for (int i = 0; i < a.GetLength(0); i++)
            {
                for (int j = 0; j < b.GetLength(1); j++)
                {
                    for (int k = 0; k < b.GetLength(0); k++)
                    {
                        r[i, j] += a[i, k] * b[k, j];
                    }
                }
            }
            return r;
        }
        public static double[,] PlusMinis(double[,] N, double[,] L, double i)
        {
            double[,] res = new double[N.GetLength(0), N.GetLength(1)];
            for (int row = 0; row < N.GetLength(0); row++)
            {
                for (int col = 0; col < N.GetLength(1); col++)
                {
                    res[row, col] = N[row, col] + L[row, col] * i;
                }
            }
            return res;
        }
        public static double[,] MatrixDuplicate(double[,] matrix)
        {
            double[,] result = new double[matrix.GetLength(0), matrix.GetLength(1)];
            for (int i = 0; i < matrix.GetLength(0); ++i)
                for (int j = 0; j < matrix.GetLength(1); ++j)
                    result[i, j] = matrix[i, j];
            return result;
        }
        public static void MatrixShow(double[,] x)
        {
            for (int i = 0; i < x.GetLength(0); i++)
            {
                for (int j = 0; j < x.GetLength(1); j++)
                {
                    Console.Write(x[i, j] + "  ");
                }
                Console.WriteLine("");
            }
            Console.WriteLine("=======================================");

        }

        public static double[,] MatrixTranspose(double[,] Matrix)
        {
            int n = Matrix.GetLength(0);
            int m = Matrix.GetLength(1);

            double[,] Result = new double[m, n];
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    Result[i, j] = Matrix[j, i];
                }
            }
            return Result;
        }
        public static double[,] MatrixInverse(double[,] matrix)
        {
            int n = matrix.GetLength(0);
            double[,] result = MatrixDuplicate(matrix);

            double[,] lum = MatrixDecompose(matrix, out int[] perm, out int toggle);
            if (lum == null)
                throw new Exception("Unable to compute inverse");

            double[] b = new double[n];
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    if (i == perm[j])
                        b[j] = 1.0;
                    else
                        b[j] = 0.0;
                }

                double[] x = HelperSolve(lum, b);

                for (int j = 0; j < n; ++j)
                    result[j, i] = x[j];
            }
            return result;
        }

        public static double[] HelperSolve(double[,] luMatrix, double[] b)
        {
            int n = luMatrix.GetLength(0);
            double[] x = new double[n];
            b.CopyTo(x, 0);

            for (int i = 1; i < n; ++i)
            {
                double sum = x[i];
                for (int j = 0; j < i; ++j)
                    sum -= luMatrix[i, j] * x[j];
                x[i] = sum;
            }

            x[n - 1] /= luMatrix[n - 1, n - 1];
            for (int i = n - 2; i >= 0; --i)
            {
                double sum = x[i];
                for (int j = i + 1; j < n; ++j)
                    sum -= luMatrix[i, j] * x[j];
                x[i] = sum / luMatrix[i, i];
            }

            return x;
        }
        public static double[,] MatrixDecompose(double[,] matrix, out int[] perm, out int toggle)
        {
            int rows = matrix.GetLength(1);
            int cols = matrix.GetLength(0);
            if (rows != cols)
                throw new Exception("Attempt to decompose a non-square m");

            int n = rows;

            double[,] result = MatrixDuplicate(matrix);

            perm = new int[n];
            for (int i = 0; i < n; ++i) { perm[i] = i; }

            toggle = 1;

            for (int j = 0; j < n - 1; ++j)
            {
                double colMax = Math.Abs(result[j, j]);
                int pRow = j;

                for (int i = j + 1; i < n; ++i)
                {
                    if (Math.Abs(result[i, j]) > colMax)
                    {
                        colMax = Math.Abs(result[i, j]);
                        pRow = i;
                    }
                }

                if (pRow != j)
                {
                    double[] rowPtr = new double[n];
                    for (int i = 0; i < n; i++) rowPtr[i] = result[pRow, i];
                    for (int i = 0; i < n; i++) result[pRow, i] = result[j, i];
                    for (int i = 0; i < n; i++) result[j, i] = rowPtr[i];

                    int tmp = perm[pRow];
                    perm[pRow] = perm[j];
                    perm[j] = tmp;

                    toggle = -toggle;
                }


                if (result[j, j] == 0)
                {
                    int goodRow = -1;
                    for (int row = j + 1; row < n; ++row)
                    {
                        if (result[row, j] != 0)
                            goodRow = row;
                    }

                    if (goodRow == -1)
                        throw new Exception("Cannot use Doolittle's method");

                    double[] rowPtr = new double[n];
                    for (int i = 0; i < n; i++) rowPtr[i] = result[goodRow, i];
                    for (int i = 0; i < n; i++) result[goodRow, i] = result[j, i];
                    for (int i = 0; i < n; i++) result[j, i] = rowPtr[i];


                    int tmp = perm[goodRow];
                    perm[goodRow] = perm[j];
                    perm[j] = tmp;

                    toggle = -toggle;
                }

                for (int i = j + 1; i < n; ++i)
                {
                    result[i, j] /= result[j, j];
                    for (int k = j + 1; k < n; ++k)
                    {
                        result[i, k] -= result[i, j] * result[j, k];
                    }
                }


            }
            return result;
        }

    }

}
