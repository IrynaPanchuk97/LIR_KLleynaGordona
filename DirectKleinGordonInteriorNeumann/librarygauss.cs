using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DirectKleinGordonInteriorNeumann
{
    public class librarygauss
    {
        public delegate double func(double x);
        public delegate double func1(double x, double z);
        public double square;
        private int No;
        private int N;

        public double IntegalTrap(double a, double b, int N, func f)
        {
            double h = (b - a) / N, s = 0, TEMP = a;
            double[] x = new double[N + 1];
            for (int i = 0; i < x.Length; i++)
            {
                x[i] = TEMP;
                if ((x[i] == a) || (x[i] == b))
                {
                    s = s + 1 / 2.0 * f(x[i]);
                }
                else
                {
                    s = s + f(x[i]);
                }
                TEMP = TEMP + h;
            }
            return (h * s);
        }

        public double[] gauss(double[,] matr, double[] b)
        {
            int N = b.Length;
            double[] x = new double[N];
            double R;
            try
            {

                for (int q = 0; q < N; q++)
                {
                    R = 1 / matr[q, q];
                    matr[q, q] = 1;
                    for (int j = q + 1; j < N; j++)
                        matr[q, j] *= R;
                    b[q] *= R;
                    for (int k = q + 1; k < N; k++)
                    {
                        R = matr[k, q];
                        matr[k, q] = 0;
                        for (int j = q + 1; j < N; j++)
                            matr[k, j] = matr[k, j] - matr[q, j] * R;
                        b[k] = b[k] - b[q] * R;
                    }
                }
            }
            catch (DivideByZeroException)
            {

                return x;
            }

            for (int q = N - 1; q >= 0; q--)
            {
                R = b[q];
                for (int j = q + 1; j < N; j++)
                    R -= matr[q, j] * x[j];
                x[q] = R;
            }
            return x;
        }

    }
}