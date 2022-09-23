#ifndef SOLVER_STUDENT_H
#define SOLVER_STUDENT_H

#include <iostream>
#include "Legendre.h"

/*!
 * Essa é uma classe comum em simuladores, e é responsável por "resolver" integrais e derivadas (mais especificamente, sistemas de equações diferenciais).
 * Seus métodos incluem, principalmente, getters e setters de parâmetros e métodos para integrar ou derivar funções com diferentes quantidades de parâmetros (fp1, fp2, ...).
 * O parâmetro p1 é o valor no qual a função está sendo avaliada, e p2, p3, ... são os parâmetros das funções.
 * Como as definições dos tipos Solver_if::f1p, Solver_if::f2p, etc não estão disponíveis (pois Solver_if não é apresentado), eles são apresentados a seguir:
 *  typedef double (*f1p)(double);
 *   typedef double (*f3p)(double, double, double);
 *   typedef double (*f4p)(double, double, double, double);
 *   typedef double (*f5p)(double, double, double, double, double);
 * , ou seja, são basicamente ponteiros para funções que retornam um double (o valor da função) e que possuem um (fp1), dois (fp2) ou mais parâmetros de entrada.
 */
class SolverStudent {
public:
    SolverStudent() {}

    /*!
    * Esse método deve setar um atributo privado com o tamanho mínimo do passo (h_min) a ser dado numa integração ou derivação numérica
    */
	virtual void setMinimumStepSize(double e) { 
        minimumStepSize = e;
	}

    /*!
    * Esse método deve retornar o atributo privado que informa o tamanho mínimo do passo (h_min) a ser dado numa integração ou derivação numérica
    */
	virtual double getMinimumStepSize() { 
	    return minimumStepSize;
	}

    /*!
    * Esse método deve setar um atributo privado com a quantidade máxima de passos a serem dados numa integração ou derivação numérica
    */
	virtual void setMaxSteps(double steps) { 
	    maxSteps = steps;
	}
	
    /*!
    * Esse método deve retornar o atributo privado que informa a quantidade máxima de passos a serem dados numa integração ou derivação numérica
    */
	virtual double getMaxSteps() { 
	    return maxSteps;
	}

    /*!
    * Todos os métodos "integrate" abaixo devem calcular a integral definida de "min" até "max" da função "f" utilizando o método da quadratura gaussina, estudado em cálculo numérico.
    * O intervalode "min" até "max" deve ser dividido na quantidade máxima de passos. 
    * Se isso resultar numa largura de passo igual ou maior que a largura mínima de um passo, tudo bem. 
    * Senão, deve-se ajustar a quantidade de passos para que a largura mínima seja a menor possível tal que não seja menor que a largura mínima.
    * Exemplo: Se o intervalo for 10, a quantidade máxima de passos por 20 e a largura mínima for 1.5, então inicialmente coloca-se 20 passos, o que gera passos de largura 10/20=0.5,
    *    mas 0.5 é menor que a largura mínima 1.5. Então pode-se calcular que 10/1.5=6.666... passos. Portanto, deve-se fazer os cálculos com 6 passos (a quantidade de passos é um inteiro), 
    *    o que leva a uma largura de passo 10/6=1.666...[OBS: para a quadratura gaussiana pode-se usar as raizes de um polinônio de Legendre de grau 5]
    * A funçao "f" pode ter 1, 2, ... , 5 parâmetros. Por isso há 5 métodos diferentes. 
    * Seus códigos serão praticamente idênticos. A única diferença deve ser na hora de invocar "f", depois demandará invocar f com 1, 2, ..., 5 parâmetros.
    */
	virtual double integrate(double min, double max, Solver_if::f1p f) {
	    double half = (max - min) / 2;
	    double middle = (max + min) / 2;
	    unsigned int degree = Legendre::maxN();
	    double result = 0;
	    for (int i = 0; i < degree; ++i) {
	        double root = 0;
	        double weight = 0;
	        Legendre::values(degree, i, root, weight);
	        result += weight * f(half * root + middle);
	    }
	    
	    return result * half;
	}

	virtual double integrate(double min, double max, Solver_if::f2p f, double p2) { 
	    double half = (max - min) / 2;
	    double middle = (max + min) / 2;
	    unsigned int degree = Legendre::maxN();
	    double result = 0;
	    for (int i = 0; i < degree; ++i) {
	        double root = 0;
	        double weight = 0;
	        Legendre::values(degree, i, root, weight);
	        result += weight * f(half * root + middle, p2);
	    }
	    
	    return result * half;
	}

	virtual double integrate(double min, double max, Solver_if::f3p f, double p2, double p3) { 
	    double half = (max - min) / 2;
	    double middle = (max + min) / 2;
	    unsigned int degree = Legendre::maxN();
	    double result = 0;
	    for (int i = 0; i < degree; ++i) {
	        double root = 0;
	        double weight = 0;
	        Legendre::values(degree, i, root, weight);
	        result += weight * f(half * root + middle, p2, p3);
	    }
	    
	    return result * half;
	}

	virtual double integrate(double min, double max, Solver_if::f4p f, double p2, double p3, double p4) { 
	    double half = (max - min) / 2;
	    double middle = (max + min) / 2;
	    unsigned int degree = Legendre::maxN();
	    double result = 0;
	    for (int i = 0; i < degree; ++i) {
	        double root = 0;
	        double weight = 0;
	        Legendre::values(degree, i, root, weight);
	        result += weight * f(half * root + middle, p2, p3, p4);
	    }
	    
	    return result * half;
	}

	virtual double integrate(double min, double max, Solver_if::f5p f, double p2, double p3, double p4, double p5) { 
	    double half = (max - min) / 2;
	    double middle = (max + min) / 2;
	    unsigned int degree = Legendre::maxN();
	    double result = 0;
	    for (int i = 0; i < degree; ++i) {
	        double root = 0;
	        double weight = 0;
	        Legendre::values(degree, i, root, weight);
	        result += weight * f(half * root + middle, p2, p3, p4, p5);
	    }
	    
	    return result * half;
	}

    /*!
    * Todos os métodos "derivate" abaixo devem resolver um sistema de equações diferenciais de primeira-ordem e valor inicial, utilizando o método de Runge-Kutta de 4a ordem a partir de um ponto inicial "min" de cada função de "std::vector<> f até um ponto final "max".
    * Como se trata de um problema de valor inicial, o valor de cada equação diferencial em "std::vector<> f" no ponto "min" é dado por "std::vector<> initValue".
    * A derivação de cada equação em "f", a partir do ponto "min" até o ponto final "max", deve se dar na quantidade máxima de passos e, em decorrência, disso, de certa largura de passo (h), desde que essa
        largura não seja menor que a largura mínima (MinimumStepSize).Se esse for o caso, a quantidade de passos deve ser ajustada para que a largura do passo seja a mais próxima da mínima (como indicado na integração)
    * É claro que o método de Runge-Kutta pode exigir o cálculo em pontos intermediários, como h/2.
    * Como se trata de um vector de equações diferenciais envolvendo possivelmente valores de outras equações desse sistema, deve-se utilizar o método de Runge-Kutta de 4a ordem para um sistema de equações.
    * Esse método é basicamente idêntico ao de uma única função, mas calcula k1, k2, k3, k4 para cada uma das funções (k1f1, k1f2, ..., k2f1, k2f2,...). Portanto, os valores intermediários também serão vetores std::vector<>k1, std::vector<>k2...
    * Essa função retorna um std::vector<> com os valores das derivadas finais de cada função "f" no ponto max.
    * Cada uma das funções em "f pode ter 1, 2, ... , 5 parâmetros. Por isso há 5 métodos diferentes. Perceba que todas as funções em "f" terão a mesma quantidade de parâmetros.
    * O parâmetro p1 de cada função é o ponto em que a função está sendo avaliada, e os demais parâmetros são parâmetros normais da função.
    * Seus códigos serão praticamente idênticos. A única diferença deve ser na hora de invocar "f", depois demandará invocar f com 1, 2, ..., 5 parâmetros.
    */
	virtual std::vector<double> derivate(double min, double max,  std::vector<double> initValue, std::vector<Solver_if::f2p> f) { 
        std::vector<double> res(initValue);
        double n = this->getMaxSteps();
        double h = (max - min) / n;
        double x0 = min;
        Solver_if::f2p f0 = f[0];
        double y0 = initValue[0];
        for (int j = 0; j < n; ++j) {
            double k1 = h * f0(x0, y0);
            double k2 = h * f0(x0 + 0.5 * h, y0 + 0.5 * k1);
            double k3 = h * f0(x0 + 0.5 * h, y0 + 0.5 * k2);
            double k4 = h * f0(x0 + h, y0 + k3);
            x0 = x0 + h;
            y0 = y0 + (k1 + k2 + k2 + k3 + k3 + k4) / 6;
        }
        
        res[0] = y0;
	    return res;
	}

	virtual std::vector<double> derivate(double min, double max,  std::vector<double> initValue, std::vector<Solver_if::f3p> f) { 
        std::vector<double> res(initValue);
        double n = this->getMaxSteps();
        double h = (max - min) / n;
        double x0, x1 = min;
        Solver_if::f3p f0 = f[0];
        double y0 = initValue[0];
        Solver_if::f3p f1 = f[1];
        double y1 = initValue[1];
        for (int j = 0; j < n; ++j) {
            double k01 = h * f0(x0, y0, y1);
            double k02 = h * f0(x0 + 0.5 * h, y0 + 0.5 * k01, y1 + 0.5 * k01);
            double k03 = h * f0(x0 + 0.5 * h, y0 + 0.5 * k02, y1 + 0.5 * k02);
            double k04 = h * f0(x0 + h, y0 + k03, y1 + k03);
            
            double k11 = h * f1(x1, y0, y1);
            double k12 = h * f1(x1 + 0.5 * h, y0 + 0.5 * k11, y1 + 0.5 * k11);
            double k13 = h * f1(x1 + 0.5 * h, y0 + 0.5 * k12, y1 + 0.5 * k12);
            double k14 = h * f1(x1 + h, y0 + k13, y1 + k13);
            
            x0 += h;
            x1 += h;
            y0 += (k01 + k02 + k02 + k03 + k03 + k04) / 6;
            y1 += (k11 + k12 + k12 + k13 + k13 + k14) / 6;
        }
        
        std::cout << "x0: " << x0 << std::endl;
        std::cout << "y0: " << y0 << std::endl;
        std::cout << "x1: " << x1 << std::endl;
        std::cout << "y1: " << y1 << std::endl;
        res[0] = y0;
        res[1] = y1;
	    return res;
	}

	virtual std::vector<double> derivate(double min, double max,  std::vector<double> initValue, std::vector<Solver_if::f4p> f) { 
	    std::vector<double> res(initValue);
        double n = this->getMaxSteps();
        double h = (max - min) / n;
        double x0, x1, x2 = min;
        Solver_if::f4p f0 = f[0];
        double y0 = initValue[0];
        Solver_if::f4p f1 = f[1];
        double y1 = initValue[1];
        Solver_if::f4p f2 = f[2];
        double y2 = initValue[2];
        for (int j = 0; j < n; ++j) {
            double k01 = h * f0(x0, y0, y1, y2);
            double k02 = h * f0(x0 + 0.5 * h, y0 + 0.5 * k01, y1 + 0.5 * k01, y2 + 0.5 * k01);
            double k03 = h * f0(x0 + 0.5 * h, y0 + 0.5 * k02, y1 + 0.5 * k02, y2 + 0.5 * k02);
            double k04 = h * f0(x0 + h, y0 + k03, y1 + k03, y2 + k03);
            
            double k11 = h * f1(x1, y0, y1, y2);
            double k12 = h * f1(x1 + 0.5 * h, y0 + 0.5 * k11, y1 + 0.5 * k11, y2 + 0.5 * k11);
            double k13 = h * f1(x1 + 0.5 * h, y0 + 0.5 * k12, y1 + 0.5 * k12, y2 + 0.5 * k12);
            double k14 = h * f1(x1 + h, y0 + k13, y1 + k13, y2 + k13);
            
            double k21 = h * f2(x2, y0, y1, y2);
            double k22 = h * f2(x2 + 0.5 * h, y0 + 0.5 * k21, y1 + 0.5 * k21, y2 + 0.5 * k21);
            double k23 = h * f2(x2 + 0.5 * h, y0 + 0.5 * k22, y1 + 0.5 * k22, y2 + 0.5 * k22);
            double k24 = h * f2(x2 + h, y0 + k23, y1 + k23, y2 + k23);
            
            x0 += h;
            x1 += h;
            x2 += h;
            y0 += (k01 + k02 + k02 + k03 + k03 + k04) / 6;
            y1 += (k11 + k12 + k12 + k13 + k13 + k14) / 6;
            y2 += (k21 + k22 + k22 + k23 + k23 + k24) / 6;
        }

        res[0] = y0;
        res[1] = y1;
        res[2] = y2;
	    return res;
	}

	virtual std::vector<double> derivate(double initPoint, double endPoint,  std::vector<double> initValue, std::vector<Solver_if::f5p> f) { 
	    std::vector<double> res(initValue);
	    return res;
	}
private:
    double minimumStepSize;
    double maxSteps;
};

#endif /* SOLVER_STUDENT_H */

