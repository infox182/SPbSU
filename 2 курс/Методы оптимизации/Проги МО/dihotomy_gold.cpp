#include <iostream>
#include <math.h>
#include <locale.h>
using namespace std;
double xd = 0; //абсцисса минимума
double xg = 0; //абсцисса минимума

double f(double x,double n ,double m) {
	return  (n/exp(x) + m * x);
};
void dihotomia(double a, double b, double d, double e,double n,double m) {
	double e1 = (b - a) / 2;
	if (e1 <= e) {
		xd = (a + b) / 2;
	}
	else {
		double x1 = (a + b - d) / 2;
		double x2 = (a + b + d) / 2;

		if (f(x1,n,m) <= f(x2,n,m)) {
			b = x2;
		}
		else {
			a = x1;
		};
		
		dihotomia(a, b, e, d, n, m);
	};
};
void golden(double a, double b, double e, double x1, double x2,double n , double m) {
	double e1 = (b - a) / 2;
	if (e1 <= e) {
		xg = (a + b) / 2;
	}
	else {
		if (f(x1,n,m) <= f(x2,n,m)) {
			b = x2;
			x2 = x1;
			x1 = a + ((3 - sqrt(5)) / 2) * (b - a);
		}
		else {
			a = x1;
			x1 = x2;
			x2 = a + ((sqrt(5) - 1) / 2) * (b - a);
		};
		golden(a, b, e, x1, x2, n, m);
	};
};
int main() {
	setlocale(0, "RUS");
	double a = -1000; // интервал поиска
	double b = 1000;    //итнервал поиска
	double e = 0.0001;   //погрешность
	double d = e * 1.5;  //шаг 

	double n, m; //коэффициенты функции
	
	cout << "¬ведите коэффициенты n и m дл€ функции : n/exp(x) + m*x : " << endl;
	cin >> n >> m;
	dihotomia(a, b, e, d, n, m);
	golden(a, b, e, a + ((3 - sqrt(5)) / 2) * (b - a), a + ((sqrt(5) - 1) / 2) * (b - a),n,m);

	cout << "dihotomia" << " " << "x(min)*=" << xd << " " << "f(min)*="<< f(xd,n,m) << endl;
	cout << "golden" << " " << "x(min)*=" << xg << " " << "f(min)*="<< f(xg,n,m) ;
	cin.get();
};