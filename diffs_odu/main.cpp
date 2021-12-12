#include <iostream>
#include <vector>
#include <algorithm>

#include "math.h"

#include <fstream>
#include <iostream>

#include <map>

double q(double x) {
 //   return exp(x) * x;
    return sin(x * x) - 3.0 / x;
}

double dq_dx(double x) {
    //return exp(x) + exp(x)*x;
    return 2.0*x * cos(x * x) + 3.0 / x / x;
}

double f(double x, double y, double k = 1) {
   // return exp(x) + y;
    return dq_dx(x) + k * (y - q(x));
}
    //-2 + 3 e^x - 2 x - x^2

//double q(double x) {
//    return -2 + 3*exp(x) - 2*x - x*x;
//}
//
//double dq_dx(double x) {
//    return 3*exp(x) - 2*(1 + x);
//}
//
//double f(double x, double y, double k = 0) {
//    return x*x+y + k * (y - q(x));
//}

//double q(double x) {
//    return sin(x*x*x / 3) / cos(x*x*x / 3);
//}

//double dq_dx(double x) {
//    return 3*exp(x) - 2*(1 + x);
//}

//double f(double x, double y, double k = 0) {
//    return x*x*(y*y + 1) + k * (y - q(x));
//}


void solve1(double a, double b, const int n, std::vector<double> xs, std::vector<double> &ys, std::vector<double> &eps, std::vector<double> rs) {
 
    

    ys.clear();
    ys.resize(n + 1);
    
    eps.clear();
    eps.resize(n + 1);
    double h = (b - a) / n;
    //std::vector<double> ys(n + 1);
    ys[0] = q(a);
    for (int i = 1; i <= n; i++) {
        ys[i] = ys[i-1] + h * f(xs[i-1], ys[i-1]);
    }
    
    for (int i = 0; i <= n; i++) {
        eps[i] = abs(rs[i] - ys[i]);
    }
    
}

void solve2(double A, double B, const int N, std::vector<double> xs, std::vector<double> &W, std::vector<double> &eps, std::vector<double> rs)
{
      double /*A,B,ALPHA,H,T,W[MAX_N],*/ K1,K2,K3,K4;
    int I;
    
    W.clear();
    W.resize(N + 1);
    
    eps.clear();
    eps.resize(N + 1);

//      A = 0.0;
//      B = 2.0;
//      N = 10;

//      cout.setf(ios::fixed,ios::floatfield);
//      cout.precision(10);

      /* STEP 1 */
      double H = (B - A) / N;
      double T = A;
      W[0] = q(xs[0]); // initial value

      /* STEP 2 --- Use order 4 RK method to get w1, w2, w3 */
      /// NOTE: The "for" loop starts with I = 1.
      for (I=1; I<=3; I++)
      {
         /* STEP 3 */
         /* Compute K1, K2 RESP. */
         K1 = H*f(T, W[I-1]);
         K2 = H*f(T + H/2.0, W[I-1] + K1/2.0);
         K3 = H*f(T + H/2.0, W[I-1] + K2/2.0);
         K4 = H*f(T + H, W[I-1] + K3);

         /* STEP 4 */
         /* COMPUTE W(I) */
         W[I] = W[I-1] + 1/6.0*(K1 + 2.0*K2 + 2.0*K3 + K4);

         /* COMPUTE T(I) */
         T = A + I * H;

         /* STEP 5 */
         //cout <<"At time "<< T <<" the solution = "<< W[I] << endl;
      }

      /* STEP 6---Use Adam-Bashforth 4-step explicit method */
      for(I = 4; I <= N; I++)
      {
          K1 = 55.0*f(T, W[I-1]) - 59.0*f(T-H, W[I-2]) + 37.0*f(T-2.0*H, W[I-3]) - 9.0*f(T-3.0*H, W[I-4]);
          W[I] = W[I-1] + H/24.0*K1;

          /* COMPUTE T(I) */
          T = A + I * H;

          /* STEP 7 */
          //cout <<"At time "<< T <<" the solution = "<< W[I] << endl;
      }

    for(I = 0; I <= N; I++) {
        eps[I] = abs(q(A + I * H) - W[I]);
    }
      
}

void solve3_a(double a, double b, const int n, std::vector<double> xs, std::vector<double> &ys, std::vector<double> &eps, std::vector<double> rs) {
    double lambda = 0.001;
    
    ys.clear();
    ys.resize(n + 1);
    
    eps.clear();
    eps.resize(n + 1);
    
    double h = (b - a) / n;
    ys[0] = q(a);
    
    
    
    for (int i = 1; i <= n; i++) {
        double y = ys[i-1] + h * f(xs[i-1], ys[i-1]);
        double y_next;
    
        size_t c = 0;
        while (true) {
            c++;
            y_next = ys[i-1] + 0.5 * h * (f(xs[i-1], ys[i-1]) + f(xs[i], y));
            
            if (abs(y_next - y) / abs(y) <= lambda || c > 5000) {
                ys[i] = y_next;
                break;
            } else {
                y = y_next;
            }
        }
    }
    
    for (int i = 0; i <= n; i++) {
        eps[i] = abs(q(xs[i]) - ys[i]);
    }
}

void solve3_b(double a, double b, const int n, std::vector<double> xs, std::vector<double> &ys, std::vector<double> &eps, std::vector<double> rs) {
    double lambda = 0.001;
    
    ys.clear();
    ys.resize(n + 1);
    
    eps.clear();
    eps.resize(n + 1);
    
    double h = (b - a) / n;
    ys[0] = q(a);
    
    double k1, k2, k3, k4;
    
    for (int i = 1; i <= 3; i++) {
        k1 = h * f(xs[i-1], ys[i-1]);
        k2 = h * f(xs[i-1] + h / 2, ys[i-1] + k1 / 2);
        k3 = h * f(xs[i-1] + h / 2, ys[i-1] + k2 / 2);
        k4 = h * f(xs[i-1] + h, ys[i-1] + k3);
        
        ys[i] = ys[i-1] + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    }
    
    for (int i = 4; i <= n; i++) {
        //явный адамса
        double y = ys[i-1] + h / 24 * (55 * f(xs[i-1], ys[i-1]) - 59 * f(xs[i-2], ys[i-2]) + 37 * f(xs[i-3], ys[i-3]) - 9 * f(xs[i-4], ys[i-4]));
        double y_next;
    
        
        size_t c = 0;
        while (true) {
            c++;
            y_next = ys[i-1] + h / 24 * (9 * f(xs[i], y) + 19 * f(xs[i-1], ys[i-1]) - 5 * f(xs[i-2], ys[i-2]) + f(xs[i-3], ys[i-3]));
            
            if (abs(y_next - y) / abs(y)< lambda || c > 5000) {
                ys[i] = y_next;
                break;
            } else {
                y = y_next;
              //  std::cout << y_next << std::endl;
            }
        }
    }
    
    for (int i = 0; i <= n; i++) {
        eps[i] = abs(q(xs[i]) - ys[i]);
    }
}


void solve4_a(double a, double b, const int n, std::vector<double> xs, std::vector<double> &ys, std::vector<double> &eps, std::vector<double> rs) {
    ys.clear();
    ys.resize(n + 1);
    
    eps.clear();
    eps.resize(n + 1);
    
    double h = (b - a) / n;
    ys[0] = q(a);
    
    double k1, k2, k3;
    
    for (int i = 1; i <= n; i++) {
        k1 = h * f(xs[i-1], ys[i-1]);
        k2 = h * f(xs[i-1] + h /2, ys[i-1] + k1 / 2);
        k3 = h * f(xs[i-1] + h, ys[i-1] - k1 + 2 * k2);
        
        ys[i] = ys[i-1] + (k1 + 4 * k2 + k3) / 6;
    }
    
    for (int i = 0; i <= n; i++) {
        eps[i] = abs(q(xs[i]) - ys[i]);
    }
}

void solve4_b(double a, double b, const int n, std::vector<double> xs, std::vector<double> &ys, std::vector<double> &eps, std::vector<double> rs) {
    ys.clear();
    ys.resize(n + 1);
    
    eps.clear();
    eps.resize(n + 1);
    
    double h = (b - a) / n;
    ys[0] = q(a);
    
    double k1, k2, k3, k4;
    
    for (int i = 1; i <= n; i++) {
        k1 = h * f(xs[i-1], ys[i-1]);
        k2 = h * f(xs[i-1] + h / 2, ys[i-1] + k1 / 2);
        k3 = h * f(xs[i-1] + h / 2, ys[i-1] + k2 / 2);
        k4 = h * f(xs[i-1] + h, ys[i-1] + k3);
        
        ys[i] = ys[i-1] + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    }
    
    for (int i = 0; i <= n; i++) {
        eps[i] = abs(q(xs[i]) - ys[i]);
    }
}


void solve5(double a, double b, const int n, std::vector<double> xs, std::vector<double> &ys, std::vector<double> &eps, std::vector<double> rs) {
    double lambda = 0.001;
    
    ys.clear();
    ys.resize(n + 1);
    
    eps.clear();
    eps.resize(n + 1);
    
    double H = (b - a) / n;
    double h = H;
    ys[0] = q(a);
    
    double k1, k2, k3, k4;
    
    for (int i = 1; i <= 3; i++) {
        k1 = h * f(xs[i-1], ys[i-1]);
        k2 = h * f(xs[i-1] + h / 2, ys[i-1] + k1 / 2);
        k3 = h * f(xs[i-1] + h / 2, ys[i-1] + k2 / 2);
        k4 = h * f(xs[i-1] + h, ys[i-1] + k3);
        
        ys[i] = ys[i-1] + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    }
    
    

    
    double a0 = 25.0 / 12;
    double a1 = -48.0 / 12;
    double a2 = 36.0 / 12;
    double a3 = -16.0 / 12;
    double a4 = 3.0 / 12;
    
  
    for (int i = 4; i <= n; i++) {
        
        // Adams explicit O(h^4)
        double y = ys[i-1] + h / 24.0 * (55 * f(xs[i-1], ys[i-1]) - 59 * f(xs[i-2], ys[i-2]) + 37 * f(xs[i-3], ys[i-3]) - 9 * f(xs[i-4], ys[i-4]));
    
        double y_next;
    
        size_t c = 0;
        while (true) {
            c++;
            y_next =  ( h * f(xs[i], y) - a1 * f(xs[i-1], ys[i-1])  - a3 * f(xs[i-3], ys[i-3]) - a4 * f(xs[i-4], ys[i-4]) - a2 * f(xs[i-2], ys[i-2]) ) / a0;
       
            
            if (abs(y_next - y) < lambda || c > 5000) {
                ys[i] = y_next;
                break;
            } else {
                y = y_next;
            }
        }
    }
    
    for (int i = 0; i <= n; i++) {
        eps[i] = abs(q(xs[i]) - ys[i]);
    }
}

//void solve5(double a, double b, const int n, std::vector<double> xs, std::vector<double> &ys, std::vector<double> &eps, std::vector<double> rs) {
//    double lambda = 0.1;
//
//    ys.clear();
//    ys.resize(n + 1);
//
//    eps.clear();
//    eps.resize(n + 1);
//
//    double h = (b - a) / n;
//    ys[0] = q(a);
//
//    double k1, k2, k3, k4;
//
//    for (int i = 1; i <= 3; i++) {
//        k1 = h * f(xs[i-1], ys[i-1]);
//        k2 = h * f(xs[i-1] + h / 2, ys[i-1] + k1 / 2);
//        k3 = h * f(xs[i-1] + h / 2, ys[i-1] + k2 / 2);
//        k4 = h * f(xs[i-1] + h, ys[i-1] + k3);
//
//        ys[i] = ys[i-1] + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
//    }
//
//
//    double a0 = 25.0 / 12;
//    double a1 = -48.0/12;
//    double a2 = 36.0/12;
//    double a3 = -16.0 / 12;
//    double a4 = 3.0/12;
//
//
//    for (int i = 4; i <= n; i++) {
//        double y = ys[i-1] + h * f(xs[i-1], ys[i-1]);
//        double y_next;
//
//        while (true) {
//            y_next = h  * f(xs[i], y) - (a1 * f(xs[i-1], ys[i-1]) + a2 * f(xs[i-2], ys[i-2]) + a3 * f(xs[i-3], ys[i-3]) + a4 * f(xs[i-4], ys[i-4]));
//            y_next /= a0;
//
//            if (abs(y_next - y) < lambda) {
//                ys[i] = y_next;
//                break;
//            } else {
//                y = y_next;
//            }
//        }
//    }
//
//    for (int i = 0; i <= n; i++) {
//        eps[i] = abs(q(xs[i]) - ys[i]);
//    }
//}

const char* methodName(int i) {
    switch (i) {
        case 1:
            return "1";
            break;
        case 2:
            return "2";
            break;
        case 3:
            return "3A";
        case 4:
            return "3B";
        case 5:
            return "4A";
        case 6:
            return "4B";
        case 7:
            return "5";
        default:
            return "";
            break;
    }
}


double methodP(int i) {
    switch (i) {
        case 1:
            return 1;
            break;
        case 2:
            return 4;
            break;
        case 3:
            return 2;
        case 4:
            return 4;
        case 5:
            return 3;
        case 6:
            return 4;
        case 7:
            return 4;
        default:
            return 1;
            break;
    }
}


int main(int argc, const char * argv[]) {
    

    const double a = 0.5;
    const double b = 3;
    
    const size_t MAX_NUM = 50000;
    
    auto methods = {&solve1, &solve2, &solve3_a, &solve3_b, &solve4_a, &solve4_b/*, &solve5*/};
//
//    const double a = 0;
//    const double b = 1;
    
    std::cout << "n" << "\t\t" << "maxEps" << "\t\t" << "ro" << std::endl;
    char filename[100];
    sprintf(filename, "method_EPS_MAX.csv");
    std::ofstream file_max_eps(filename, std::ofstream::out);
    file_max_eps << std::fixed;
    file_max_eps.precision(15);
    
    file_max_eps << "N, 1, p*, r, 2, p*, r, 3A, p*, r, 3B, p*, r, 4A, p*, r, 4B, p*, r" << std::endl;
    //file_max_eps << "N, 1, p*, r, 2, p*, r, 3A, p*, r, 3B, p*, r, 4A, p*, r, 4B, p*, r, 5, p*, r" << std::endl;
    
    std::map<std::string, double> last_eps_max;
    last_eps_max["1"] = -1;
    last_eps_max["2"] = -1;
    last_eps_max["3A"] = -1;
    last_eps_max["3B"] = -1;
    last_eps_max["4A"] = -1;
    last_eps_max["4B"] = -1;
   // last_eps_max["5"] = -1;
    
    
    for (int n = 5; n < MAX_NUM; n *= 2) {
        std:: cout << std::endl << std::setw(6) << n ;
        file_max_eps << n << ", ";
        
        std::vector<double> xs(n + 1);
        
        double h = (b - a) / n;
        for (int i = 0; i <= n; i++) {
            xs[i] = a + i * h;
        }
        
        // ================================================================
        std::vector<double> ys, eps;
        
        std::vector<double> rs(n + 1);
        
      
        for (int i = 0; i <= n; i++) {
            rs[i] = q(xs[i]);
            
        }
        
        
        
        
        // ================================================================
        int count = 0;
        for (auto method : methods) {
            count++;
            method(a, b, n, xs, ys, eps, rs);
            
           
            sprintf(filename, "method_%s_%d.txt", methodName(count), n);
            std::ofstream file(filename, std::ofstream::out);
            file.precision(16);
            
            for (int i = 0; i <= n; i++) {
                file << std::fixed << xs[i] << '\t' << std::fixed << ys[i] << '\n';
            }
           
            
            double maxEps = *std::max_element(std::begin(eps), std::end(eps));
            file_max_eps << maxEps;
            
            file_max_eps << ", ";
            
            
            double prev_eps_max = last_eps_max[methodName(count)];
            last_eps_max[methodName(count)] = maxEps;
            if (prev_eps_max < 0) {
                file_max_eps << "-42, -42";
            } else {
                file_max_eps <<log2(prev_eps_max / maxEps) << ", " << abs(prev_eps_max - maxEps) / (pow(2.0, methodP(count)) - 1.0);
            }
            file_max_eps << ", ";
            
            std::cout << std::setw(16) << (abs(prev_eps_max - maxEps) < 0.001 ? h : -1) ;
            
            //std::cout << n << "\t\t" << maxEps << "\t\t" << maxEps / h/h/h/h << std::endl;
        }
    
        if (n * 2 < MAX_NUM) {
            file_max_eps << "\n";
        }
            
        
            
        
    }
    return 0;
}
