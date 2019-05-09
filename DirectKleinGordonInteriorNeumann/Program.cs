using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DirectKleinGordonInteriorNeumann
{
    class Program
    {
        static double[] yStar = { -2, 4 };
        static double alpha = Math.Exp(-10);
        const int k = 2;
        static double x1(double t) => 3*Math.Cos(t);
        static double x2(double t)=> Math.Sin(t)+ Math.Cos(t);
        static double x1d(double t)=> -3*Math.Sin(t);
        static double x2d(double t)=> Math.Cos(t) - Math.Sin(t);      
        static double x1dd(double t) => -3*Math.Cos(t);   
        static double x2dd(double t)=> -Math.Sin(t)-Math.Cos(t);      
        static double diffModule(double t, double tau)=> Math.Sqrt((x1(t) - x1(tau)) * (x1(t) - x1(tau)) + (x2(t) - x2(tau)) * (x2(t) - x2(tau)));
        static double derModule(double t)=> Math.Sqrt(x1d(t) * x1d(t) + x2d(t) * x2d(t));

        static double h(double t, double tau)
        {
            double numerator = ((x1(tau) - x1(t)) * x2d(tau) - (x2(tau) - x2(t)) * x1d(tau));
            double denominator = diffModule(t, tau);
            return numerator / denominator;
        }
        static double L(double t, double tau)=> -2 * k * alglib.besselk1(k * diffModule(t, tau)) * h(t, tau);
  
        static double L1(double t, double tau)
        {
            if (t != tau) return -k * alglib.besseli1(k * diffModule(t, tau)) * h(t, tau);
            return 0;
        }

        static double L2(double t, double tau)
        {
            if (t != tau)
            {
                return L(t, tau) - L1(t, tau) * Math.Log(4 * Math.Sin((t - tau) / 2) * Math.Sin((t - tau) / 2));
            }
            else
            {
                double numerator = -(x1d(t) * x2dd(t) - x2d(t) * x1dd(t));
                double denominator = derModule(t) * derModule(t);
                return numerator / denominator;
            }
        }
        static double D(double t, double tau)=>  2*alglib.besselk0(k * diffModule(t, tau)) * derModule(tau);
        static double D1(double t, double tau)
        {
            if (t != tau)
            {
                return -alglib.besseli0(k * diffModule(t, tau)) * derModule(tau) ;
            }
            return -derModule(t);
        }

        static double D2(double t, double tau)
        {
            if (t != tau)
            {
                return D(t, tau) - D1(t, tau) * Math.Log(4 * Math.Sin((t - tau) / 2) * Math.Sin((t - tau) / 2));
            }
            return 2*derModule(t) * (Math.Log(2 / (k * derModule(t)))  - EulerConst);
        }
        static double R(int j, double t, int n, double[] ti)
        {
            double sum = 0;
            for (int m = 1; m <= n - 1; m++)
            {
                sum = sum + Math.Cos(m * (t - ti[j])) / m;
            }
            sum = sum + Math.Cos(n * (t - ti[j])) / (2 * n);
            return -sum / n;
        }
        static double exact(double[] x)=> alglib.besselk0(k * module(x, yStar)) / (2 * Math.PI);
        static double module(double[] x, double[] y)=> Math.Sqrt((x[0] - y[0]) * (x[0] - y[0]) + (x[1] - y[1]) * (x[1] - y[1]));
        static double g(double t) => -k * alglib.besselk1(k * module(new double[] { x1(t), x2(t) }, yStar)) * ((x1(t) - yStar[0]) * x2d(t) - (x2(t) - yStar[1]) * x1d(t)) / (2 * Math.PI * derModule(t) * module(new double[] { x1(t), x2(t) }, yStar));
        

        static double G(int i,double [] ti,int n)
        {
            double sum = 0;
            for (int j=0;j<2*n;j++)
            {
                sum = sum + g(ti[j]) * (D1(ti[i], ti[j]) * R(j, ti[i], n, ti) + D2(ti[i], ti[j]) / (2 * n));
            }
            return sum;
        }
        static double tildaKernel(double[] x, double tau)
        {
            double numerator = -k*alglib.besselk1(k* module(x, new double[] { x1(tau), x2(tau) })) *((x1(tau) - x[0]) * x2d(tau) - (x2(tau) - x[1]) * x1d(tau));
            double denominator = module(x, new double[] { x1(tau), x2(tau) }) ;
            return numerator / denominator;
        }
        static double numericalSolution(double[] x, int n, double[] density, double[] ti)
        {
            double Un = 0;
            for (int i = 0; i < 2 * n; i++)
            {
                Un = Un + g(ti[i]) *alglib.besselk0(k*module(x, new double[] { x1(ti[i]), x2(ti[i]) })) * derModule(ti[i]);
            }
            for (int i = 0; i < 2 * n; i++)
            {
                Un = Un - density[i] * tildaKernel(x, ti[i]);
            }
            Un = Un / (2 * n);
            return Un;
        }
        const double EulerConst = 0.5772156649;
        
        static void Main(string[] args)
        {
            double[] testPoint = { 0,0};
            int n = 2;
            for (int k = 1; k <= 7; k++)
            {
                double[] ti = new double[2 * n];
                for (int i = 0; i < 2 * n; i++) ti[i] = Math.PI * i / n;
                

                double[,] Matrix = new double[2 * n, 2 * n];
                double[,] Matrix2 = new double[2 * n, 2 * n];
                for (int i = 0; i < 2 * n; i++)
                {
                    for (int j = 0; j < 2 * n; j++)
                    {
                        Matrix[i, j] = R(j, ti[i], n, ti) * L1(ti[i], ti[j]) + L2(ti[i], ti[j]) / (2 * n);
                        if (i == j) Matrix[i, j] = Matrix[i, j] + 1;
                    }
                }
                double[] rightPart = new double[2 * n];
                for (int i = 0; i < 2 * n; i++)
                {
                    rightPart[i] = G(i,ti,n);
                }
                rightPart = MatrixHelper.MatrixVector2(MatrixHelper.MatrixTranspose(Matrix), rightPart);

                Matrix = MatrixHelper.Multiplication(MatrixHelper.MatrixTranspose(Matrix), Matrix);
                Matrix = MatrixHelper.PlusMinis(Matrix, MatrixHelper.MatrixIdentity(2 * n, alpha), 1);
                librarygauss executor = new librarygauss();
                var density = executor.gauss(Matrix, rightPart);
                Console.WriteLine("n = " + n + "    error = " + Math.Abs(numericalSolution(testPoint, n, density, ti) - exact(testPoint)));
                n = n * 2;
                
            }
            Console.ReadKey();
        }
    }
}
