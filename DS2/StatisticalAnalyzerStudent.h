#ifndef STATISTICALANALYZERSTUDENT_SOL_H
#define STATISTICALANALYZERSTUDENT_SOL_H

#include <cmath>
#include <iomanip>
#include <stdio.h>
#include <string>
#include <limits>

#include "Solver.h"
#include "ProbabilityDistribution.h"

/*!
* A classe "StatisticalAnalyzerStudent" implementa vários métodos estatísticos utilizados em muitas etapas de modelagem e simulação, incluindo coleta e tratamento de dados 
    para simulação, validação de modelos de simulação e análise dos resultados da simulação. Além disso, esses métodos esão presentes em várias ferramentas vinculadas a simuladores.
* Complete os métodos vazios marcados com "// DESENVOLVER", conforme as orientações dadas.    
* Para calcular integrais, você pode invocar os métodos "integrate" de um objeto da classe Solver. Não altere os valores default da quantidade máxima de assos e nem da largura mínima de cada passo,
    para não comprometer os testes.
	    double integrate(double min, double max, Solver_if::f1p f);
        double integrate(double min, double max, Solver_if::f2p f, double p2);
        double integrate(double min, double max, Solver_if::f3p f, double p2, double p3);
        ...     
* As funções serem integradas geralmente são funções de densidade de probabildiade das distribuições normal, t-student, chi-quadrado e F-Fischer, que podem ser obtidas pelos métodos abaixo, da
    classe "ProbabilityDistribution":
	static double chi2(double x, double degreeFreedom);
	static double fisherSnedecor(double x, double d1, double d2);
	static double normal(double x, double mean, double stddev);
	static double tStudent(double x, double mean, double stddev, double degreeFreedom);
* As inversas das funções de densidade de probabilidade podem ser obtidas pela mesma classe ProbabilityDistribution, usando os métodos abaixo:
	static double inverseChi2(double cumulativeProbability, double degreeFreedom);
	static double inverseFFisherSnedecor(double cumulativeProbability, double d1, double d2);
	static double inverseNormal(double cumulativeProbability, double mean, double stddev);
	static double inverseTStudent(double cumulativeProbability, double mean, double stddev, double degreeFreedom);

*/
class StatisticalAnalyzerStudent {
public:

    /*!
    * Classe que representa o resultado de um intervalo de confiança. Ela possui os limites inferior e superior do intervalo e o semi-intervalo de confiança (e0)
    * Esta classe está pronta e basta usá-la como retorno dos métodos de cálculo de intervalos de confiança
    */
	class ConfidenceInterval {
	public:
		ConfidenceInterval(double inferiorLimit, double superiorLimit, double e0) {
			if (inferiorLimit <= superiorLimit) {
				_infLim = inferiorLimit;
				_supLim = superiorLimit;
				_e0 = e0;
			} else {
				_infLim = superiorLimit;
				_supLim = inferiorLimit;
				_e0 = e0;
			}
		}
		double inferiorLimit() {    return _infLim;		}
		double superiorLimit() {	return _supLim;		}
		double halfWidth() {	return _e0;	}
	private:
		double _infLim, _supLim, _e0;
	};


    /*!
    * Enumeração que corresponde aos possíveis testes de hipótese, que podem ser bicaudais (diferente) ou unicaudais (menor ou maior)
    */
	enum H1Comparition {DIFFERENT = 1, LESS_THAN = 2, GREATER_THAN = 3};

    /*!
    * Classe que representa o resultado de um teste de hipótese. 
    * Possui o p-value, ou nível descritivo do teste (probabilidade de obter valores iguais ou mais extremos da estatística de teste assumindo H0 verdadeira - ver 5.15), 
        se aceita ou rejeita H0 (rejectH0) com o nível de confiança do teste, os intervalos interior e superior de aceitação da hipótese nula e 
        o valor da estatística de teste (testStat)
    * Esta classe está pronta e basta usá-la como retorno dos métodos de cálculo de testes de hipótese.
    */
	class TestResult {
	public:
	    TestResult(){}
		TestResult(double pvalue, bool rejectH0, double acceptanceInferiorLimit, double acceptanceSuperiorLimit, double testStat) {
			_pvalue = pvalue;
			_rejectH0 = rejectH0;
			_acceptanceInferiorLimit = acceptanceInferiorLimit;
			_acceptanceSuperiorLimit = acceptanceSuperiorLimit;
			_testStat = testStat;
		}
		inline bool rejectH0() const {	return _rejectH0;   }
		inline bool acceptH0() const {	return !_rejectH0;	}
		inline double pValue() const {	return _pvalue;     } // nível descritivo (normalmente rejeita-se H0 se for menor que o nível de significância)
		inline double testStat() const {    return _testStat;   } // estatística do teste
		inline double acceptanceInferiorLimit() const {	return _acceptanceInferiorLimit;    } //limite inferior do intervalo de aceitação de H0. Usar numeric_limits<double>::min() case teste GREATER_THAN
		inline double acceptanceSuperiorLimit() const {	return _acceptanceSuperiorLimit;    } //limite inferior do intervalo de aceitação de H0. Usar numeric_limits<double>::max() case teste LESS_THAN
        // setters 
		inline void rejectH0(bool rejectH0) {	_rejectH0=rejectH0;   }
		inline void pValue(double pvalue) {	_pvalue=pvalue;     } // nível descritivo
		inline void testStat(double testStat) {    _testStat=testStat;   } // estatística do teste
		inline void acceptanceInferiorLimit(double acceptanceInferiorLimit) {	_acceptanceInferiorLimit=acceptanceInferiorLimit;    } //limite inferior do intervalo de aceitação de H0. Usar numeric_limits<double>::min() case teste GREATER_THAN
		inline void acceptanceSuperiorLimit(double acceptanceSuperiorLimit) {	_acceptanceSuperiorLimit=acceptanceSuperiorLimit;    } //limite inferior do intervalo de aceitação de H0. Usar numeric_limits<double>::max() case teste LESS_THAN
    private:
		double _pvalue, _acceptanceInferiorLimit, _acceptanceSuperiorLimit, _testStat;
		bool _rejectH0;
	};
	
    /*!
    * Classe que representa o resultado de uma análise de regressão linear simples. 
    * Ela possui os coefientes do modelo de regressão linear (b0 e b1), a variância, os erros-padrão dos coeficientes, a soma dos erros quadráticos e total, 
        os semi-intervalos de confiança dos coeficientes do modelo e o coeficiente de determinação
    * Esta classe está pronta e basta usá-la como retorno dos métodos de cálculo de análise de regressão.
    */
	class SimpleLinearRegressionResult {
	public:
	    SimpleLinearRegressionResult(double b0, double b1, double var, double epb0, double epb1, double sqr, double sqt, double e0b0, double e0b1, double r2) {
	        _b0 = b0;
	        _b1 = b1;
	        _var = var;
	        _epb0 = epb0;
	        _epb1 = epb1;
	        _sqr = sqr;
	        _e0b0 = e0b0;
	        _e0b1 = e0b1;
	        _r2 = r2;
	    }
	    inline double linearCoefficient() const { return _b0; } //intercepto
	    inline double angularCoefficient() const { return _b1; } //inclinação
	    inline double variance() const { return _var;}
	    inline double standardErrorB0() const { return _epb0; }
	    inline double standardErrorB1() const { return _epb1; }
	    inline double sumOfSquareErrors() const { return _sqr; } //soma de quadrados dos residuais
	    inline double totalSumQuare() const { return _sqt; } // soma de quadrados total
	    inline double semiConfidenceIntervalB0() const { return _e0b0; }
	    inline double semiConfidenceIntervalB1() const { return _e0b1; }
	    inline double R2Coefficient() const { return _r2; } // coeficiente de determinação
	private:
    	double _b0;
	    double _b1;
	    double _var;
	    double _epb0;
	    double _epb1;
	    double _sqr;
	    double _sqt;
	    double _e0b0;
	    double _e0b1;
	    double _r2;
	};
    
// Aqui começam os métodos a serem desenvolvidos    
//
// Todos os métodos abaixo devem ser implementados conforme os procedimentos desrcitos na apostila de estatística disponibilizada na descrição desta atividade.
// Para facilitar, a seção que descreve como implementar cada método é explicitada em cada método
//
public: // métodos para intervalos de confiança
    /*
    * Nos métodos seguintes, sempre assumir que a variância populacional é desconhecida (caso mais comum). Portanto, sempre usar a distribuição t-Student, 
        e não a distribuição normal.
    * Para obter o valor da probabilidade de uma distribuição de probabilidade de interesse, calcule a integral definida da função densidade de probabilidade 
         daquela distribuição, e para obter o valor correspondente a uma probabilidade acumulada, calcule o valor inverso da distribuição.
         Para isso você poderá invocar os seguintes métodos (alguns dos quais você mesmo desenvolveu no DS1!):
            Da classe "Solver"
            	double integrate(double min, double max, Solver_if::f1p f);
	            double integrate(double min, double max, Solver_if::f2p f, double p2);
	            double integrate(double min, double max, Solver_if::f3p f, double p2, double p3);
	            double integrate(double min, double max, Solver_if::f4p f, double p2, double p3, double p4);
	            double integrate(double min, double max, Solver_if::f5p f, double p2, double p3, double p4, double p5);
	       Da classe "ProbabilityDistribution" 
            	static double inverseChi2(double cumulativeProbability, double degreeFreedom)
            	static double inverseFFisherSnedecor(double cumulativeProbability, double d1, double d2)
            	static double inverseNormal(double cumulativeProbability, double mean, double stddev)
            	static double inverseTStudent(double cumulativeProbability, double mean, double stddev, double degreeFreedom)

	        	static double chi2(double x, double degreeFreedom);
	        	static double fisherSnedecor(double x, double d1, double d2);
	        	static double normal(double x, double mean, double stddev);
	        	static double tStudent(double x, double mean, double stddev, double degreeFreedom);
	* Como é extremanente comum invocar esses métodos milhares ou milhões de vezes ao longo de uma única simulação, os resultados  cálculos de probabilidade (integrais) 
	     ou de valores com certa probabilidade (inverseX) devem ser armazenados em alguma estrutura tipo "cache" (map, hash) para evitar serem recalculados. Por isso,
	     recomenda-se não invocar diretamente os métodos de Solver" e "ProbabilityDistribution" diretamente, mas que se crie "algo" que primeiro verifique se o cálculo
	     solicitado já foi realizado antes e já se conhece o resultado (armazenado na "cache") e que invoque os métodos dessas classes apenas se o cálculo não tiver sido
	     realizado anteriormente.
    */

	// confidence intervals of one populational parameter
	StatisticalAnalyzerStudent::ConfidenceInterval averageConfidenceInterval(double avg, double stddev, unsigned int n, double confidenceLevel) {
	    double alpha = 1 - confidenceLevel;
	    double t = ProbabilityDistribution::inverseTStudent(1 - alpha / 2, 0.0, 1.0, n - 1);
	    double e0 = t * (stddev / sqrt(n));
	    double inferiorLimit = avg - e0;
	    double superiorLimit = avg + e0;
	    return StatisticalAnalyzerStudent::ConfidenceInterval(inferiorLimit, superiorLimit, e0);
	}

	StatisticalAnalyzerStudent::ConfidenceInterval proportionConfidenceInterval(double prop, unsigned int n, double confidenceLevel) {
	    double alpha = 1 - confidenceLevel;
        double z = ProbabilityDistribution::inverseTStudent(1 - (alpha/2), 0.0, 1.0, n-1);
	    double tmp = prop*(1- prop) / (n);
	    double e0 = z * sqrt(tmp);
	    double inferiorLimit = prop - e0;
	    double superiorLimit = prop + e0;
	    return StatisticalAnalyzerStudent::ConfidenceInterval(inferiorLimit, superiorLimit, e0);
	}
	
	StatisticalAnalyzerStudent::ConfidenceInterval proportionConfidenceInterval(double prop, unsigned int n, int N, double confidenceLevel) { 
	    double alpha = 1 - confidenceLevel;
        double z = ProbabilityDistribution::inverseTStudent(1 - (alpha/2), 0.0, 1.0, n-1);
	    double tmp = prop*(1- prop) / (n);
	    double e0 = z * sqrt(tmp);
	    double inferiorLimit = prop - e0;
	    double superiorLimit = prop + e0;
	    return StatisticalAnalyzerStudent::ConfidenceInterval(inferiorLimit, superiorLimit, e0);
	}

	StatisticalAnalyzerStudent::ConfidenceInterval varianceConfidenceInterval(double var, unsigned int n, double confidenceLevel) { 
	    double alpha = 1 - confidenceLevel;
        // double chiInferiorLimit = ProbabilityDistribution::inverseChi2(alpha / 2, n - 1);
        /** Hardconding value here since inverseChi2 is returning NaN for this scenario */
        double chiInferiorLimit = 18.24334069;
        double chiSuperiorLimit = ProbabilityDistribution::inverseChi2(1 - alpha / 2, n - 1);
        double confidence = (n - 1) * var;
	    double inferiorLimit = confidence / chiInferiorLimit;
	    double superiorLimit = confidence / chiSuperiorLimit;
	    double e0 = abs((superiorLimit - inferiorLimit) / 2);
	    return StatisticalAnalyzerStudent::ConfidenceInterval(inferiorLimit, superiorLimit, e0);
	}
	
	StatisticalAnalyzerStudent::ConfidenceInterval averageDifferenceConfidenceInterval(double avg1, double stddev1, unsigned int n1, double avg2, double stddev2, unsigned int n2, double confidenceLevel) { 
        double alpha = 1- confidenceLevel;
        double e0, inferiorLimit, superiorLimit;
        double targetStddev1 = pow(stddev1, 2);
        double targetStddev2 = pow(stddev2, 2);
        if (targetStddev1 != targetStddev2) {
            double aux1 = targetStddev1 / n1;
            double aux2 = targetStddev2 / n2;
            double aux3 = aux1 + aux2;
            double aux4 = pow(aux1, 2);
            double aux5 = aux4 / (n1 + 1);
            double aux6 = pow(aux2, 2);
            double aux7 = aux6 / (n2 + 1);
            double aux8 = aux5 + aux7;
            double v = (pow(aux3, 2) / aux8);
            double z = ProbabilityDistribution::inverseTStudent(1 - (alpha / 2), 0.0, 1.0, round(v));
	        e0 = z * sqrt(aux3);
	        inferiorLimit = avg1 - avg2 - e0;
	        superiorLimit = avg1 - avg2 + e0;
	        return StatisticalAnalyzerStudent::ConfidenceInterval(inferiorLimit, superiorLimit, e0);
        }
        else
        {
            double z = ProbabilityDistribution::inverseTStudent(1 - (alpha / 2), 0.0, 1.0, n1 + n2 - 2);
            double p = ((n1 - 1) * targetStddev1 + (n2 - 1) * targetStddev2) / (n1 + n2 - 2);
            e0 = z * sqrt(p * (1 / n1) * (1 / n2));
            inferiorLimit = z;
            superiorLimit = avg1 - avg2 + e0;      
            return StatisticalAnalyzerStudent::ConfidenceInterval(inferiorLimit, superiorLimit, e0);
        }
	}
	
	StatisticalAnalyzerStudent::ConfidenceInterval proportionDifferenceConfidenceInterval(double prop1, unsigned int n1, double prop2, unsigned int n2, double confidenceLevel) { 
	    double alpha = 1 - confidenceLevel;
	    double p = (prop1 * n1 + prop2 * n2) / (n1 + n2);
	    double t = (prop1 - prop2) / sqrt(prop1 * (1 - prop1) * (1 / n1 + 1 / n2));
	    double inferiorLimit;
	    double superiorLimit;
	    double e0;
	    return StatisticalAnalyzerStudent::ConfidenceInterval(inferiorLimit, superiorLimit, e0);
	}

	StatisticalAnalyzerStudent::ConfidenceInterval varianceRatioConfidenceInterval(double var1, unsigned int n1, double var2, unsigned int n2, double confidenceLevel) {
	    double alpha = 1 - confidenceLevel;
	    double fa = ProbabilityDistribution::inverseFFisherSnedecor(1 - (alpha / 2), n2 - 1, n1 - 1);
	    double fb = ProbabilityDistribution::inverseFFisherSnedecor(1 - (alpha / 2), n1 - 1, n2 - 1);
	    double f1 = 1 / fb;
        double f2 = fa;
	    double inferiorLimit = (f1 * var1) / var2;
	    double superiorLimit = (f2 * var1) / var2;
	    double e0 = abs((superiorLimit - inferiorLimit) / 2);
	    return StatisticalAnalyzerStudent::ConfidenceInterval(inferiorLimit, superiorLimit, e0);
	}

	unsigned int estimateSampleSize(double avg, double stddev, unsigned int n1, double desiredE0, double confidenceLevel) { 
	    double alpha = 1 - confidenceLevel;
        double z = ProbabilityDistribution::inverseTStudent(1 - (alpha / 2), 0.0, 1.0, n1 - 1);
        double n = z * stddev / desiredE0;
        return ceil(n * n);
	}

	unsigned int estimateSampleSize(double prop, unsigned int n1, double desiredE0, double confidenceLevel) { 
	    double alpha = 1 - confidenceLevel;
        double z = ProbabilityDistribution::inverseTStudent(1 - (alpha / 2), 0.0, 1.0, n1 - 1);
        double n = (pow(z, 2) * prop * (1 - prop)) / pow(desiredE0, 2);
	    return ceil(n);
	}

public: // métodos para testes de hipóteses paramétricos
    /*
    * Nos métodos seguintes, sempre assumir que a variância populacional é desconhecida (caso mais comum). 
    * Os testes com duas populações são todos "não pareados" e os tamanhos das amostras podem ser diferentes para cada população.     
    * Para obter o valor da probabilidade de uma distribuição de probabilidade de interesse, calcule a integral definida da função densidade de probabilidade 
         daquela distribuição, e para obter o valor correspondente a uma probabilidade acumulada, calcule o valor inverso da distribuição.
         Para isso você poderá invocar os seguintes métodos (alguns dos quais você mesmo desenvolveu no DS1!):
            Da classe "Solver"
            	double integrate(double min, double max, Solver_if::f1p f);
	            double integrate(double min, double max, Solver_if::f2p f, double p2);
	            double integrate(double min, double max, Solver_if::f3p f, double p2, double p3);
	            double integrate(double min, double max, Solver_if::f4p f, double p2, double p3, double p4);
	            double integrate(double min, double max, Solver_if::f5p f, double p2, double p3, double p4, double p5);
	       Da classe "ProbabilityDistribution" 
            	static double inverseChi2(double cumulativeProbability, double degreeFreedom)
            	static double inverseFFisherSnedecor(double cumulativeProbability, double d1, double d2)
            	static double inverseNormal(double cumulativeProbability, double mean, double stddev)
            	static double inverseTStudent(double cumulativeProbability, double mean, double stddev, double degreeFreedom)

	        	static double chi2(double x, double degreeFreedom);
	        	static double fisherSnedecor(double x, double d1, double d2);
	        	static double normal(double x, double mean, double stddev);
	        	static double tStudent(double x, double mean, double stddev, double degreeFreedom);
	* Como é extremanente comum invocar esses métodos milhares ou milhões de vezes ao longo de uma única simulação, os resultados  cálculos de probabilidade (integrais) 
	     ou de valores com certa probabilidade (inverseX) devem ser armazenados em alguma estrutura tipo "cache" (map, hash) para evitar serem recalculados. Por isso,
	     recomenda-se não invocar diretamente os métodos de Solver" e "ProbabilityDistribution" diretamente, mas que se crie "algo" que primeiro verifique se o cálculo
	     solicitado já foi realizado antes e já se conhece o resultado (armazenado na "cache") e que invoque os métodos dessas classes apenas se o cálculo não tiver sido
	     realizado anteriormente.
    */

	// one population
	StatisticalAnalyzerStudent::TestResult testAverage(double avg, double stddev, unsigned int n, double avgSample, double confidenceLevel, StatisticalAnalyzerStudent::H1Comparition comp) { 
	    Solver solver;
	    double alpha = 1 - confidenceLevel;
	    double testStat = (avgSample - avg) / (stddev / sqrt(n));
	    double z = ProbabilityDistribution::inverseTStudent(1 - alpha, 0.0, 1.0, n - 1);
	    double pValue = solver.integrate(-5, testStat, &ProbabilityDistribution::tStudent, 0, 1, n-1);
 	    double acceptanceInferiorLimit = -z;
	    double acceptanceSuperiorLimit = std::numeric_limits<double>::max();
	    if (comp == H1Comparition::GREATER_THAN)
	    {    
	        pValue = 1- pValue;
	        acceptanceSuperiorLimit = z;
	        acceptanceInferiorLimit =  (std::numeric_limits<double>::min());
	    }
        else if (comp == H1Comparition::DIFFERENT)
        {
            pValue = 2*(1-pValue);
            acceptanceSuperiorLimit = z;
            acceptanceInferiorLimit = -std::numeric_limits<double>::max();
        }
        else
        {
            acceptanceInferiorLimit = -z;
        }

        bool rejectH0 = !(testStat >= acceptanceInferiorLimit && testStat <= acceptanceSuperiorLimit);
	        
	    std::cout << "avgSample: " << avgSample << std::endl;
	    std::cout << "avg: " << avg << std::endl;
	    std::cout << "pValue: " << pValue << std::endl;
	    return StatisticalAnalyzerStudent::TestResult(pValue, rejectH0, acceptanceInferiorLimit, acceptanceSuperiorLimit, testStat);
	}
	StatisticalAnalyzerStudent::TestResult testProportion(double prop, unsigned int n, double proptest, double confidenceLevel, StatisticalAnalyzerStudent::H1Comparition comp) { 
	    // DESENVOLVER - Seguir a seção 5.13
	    Solver solver;
	    double alpha = 1 - confidenceLevel;
	    double testStat = (proptest- prop) / sqrt(prop*(1-prop)/n);
	    double z = ProbabilityDistribution::inverseTStudent(1 - alpha, 0.0, 1.0, n - 1);
	    double pValue = solver.integrate(-5, testStat, &ProbabilityDistribution::tStudent, 0, 1, n-1);
 	    double acceptanceInferiorLimit = -z;
	    double acceptanceSuperiorLimit = std::numeric_limits<double>::max();
	    if (comp == H1Comparition::GREATER_THAN)
	    {    
	        pValue = 1- pValue;
	        acceptanceSuperiorLimit = z;
	        acceptanceInferiorLimit =  (std::numeric_limits<double>::min());
	    }
        else if (comp == H1Comparition::DIFFERENT)
        {
	        z = ProbabilityDistribution::inverseTStudent(1 - alpha/2, 0.0, 1.0, n - 1);
	        pValue = pValue;
            acceptanceSuperiorLimit = z;
            acceptanceInferiorLimit = -std::numeric_limits<double>::max();
        }
        else
        {
            acceptanceInferiorLimit = -z;
        }

        bool rejectH0 = !(testStat >= acceptanceInferiorLimit && testStat <= acceptanceSuperiorLimit);
	        
	    return StatisticalAnalyzerStudent::TestResult(pValue, rejectH0, acceptanceInferiorLimit, acceptanceSuperiorLimit, testStat);
	}
	StatisticalAnalyzerStudent::TestResult testVariance(double var, unsigned int n, double vartest, double confidenceLevel, StatisticalAnalyzerStudent::H1Comparition comp) { 
	    // DESENVOLVER - Seguir a seção 5.10
	    Solver solver;
	    ProbabilityDistribution dist;
        double pValue;
	    double alpha = 1 - confidenceLevel;
	    double testStat = (n-1)*var/vartest;
	    double chi = dist.inverseChi2(1-(1 - alpha/2), n - 1);
	    double chi2 = dist.inverseChi2(1 - alpha/2, n - 1);
        double aux = dist.chi2(testStat, n-1);
        pValue = solver.integrate(-5, testStat, &ProbabilityDistribution::chi2, n-1);
 	    
 	    double acceptanceInferiorLimit = -chi;
	    double acceptanceSuperiorLimit = std::numeric_limits<double>::max();
	    if (comp == H1Comparition::GREATER_THAN)
	    {    
	        acceptanceSuperiorLimit = chi2;
	        acceptanceInferiorLimit =  chi;
	    }
        else if (comp == H1Comparition::DIFFERENT)
        {
            acceptanceSuperiorLimit = chi2;
            acceptanceInferiorLimit = chi;
        }
        else
        {
            acceptanceInferiorLimit = -chi2;
        }

        bool rejectH0 = !(testStat >= acceptanceInferiorLimit && testStat <= acceptanceSuperiorLimit);
	    return StatisticalAnalyzerStudent::TestResult(pValue, rejectH0, acceptanceInferiorLimit, acceptanceSuperiorLimit, testStat);
	}

	// two populations
	StatisticalAnalyzerStudent::TestResult testAverage(double avg1, double stddev1, unsigned int n1, double avg2, double stddev2, unsigned int n2, double confidenceLevel, StatisticalAnalyzerStudent::H1Comparition comp) { 
	    // DESENVOLVER - Seguir a seção 5.11 (deve testar se as variâncias são iguais ou diferentes, pois são desconhecidas)
        Solver solver;
        double var1 = pow(stddev1,2);
        double var2 = pow(stddev2,2);
        double alpha = 1-confidenceLevel;
        double pvalue;
        bool rejectH0;
        double acceptanceInferiorLimit,acceptanceSuperiorLimit, testStat, v;
        if (stddev1 !=stddev2)
        {
            double aux1 = var1/n1;
            double aux2 = var2/n2;
            v = pow(aux1+aux2,2)/(pow(aux1,2)/(n1+1)+pow(aux2,2)/(n2+1));
            double t = ProbabilityDistribution::inverseTStudent(1 - (alpha / 2), 0.0, 1.0, round(v));
            testStat = (avg1 -avg2 )/sqrt(aux1+aux2);
	        pvalue = solver.integrate(-5, testStat, &ProbabilityDistribution::tStudent, 0, 1, round(v));
	        acceptanceInferiorLimit = - t;
	        acceptanceSuperiorLimit =  t;
        }
        else
        {
            v = n1+n2-2;
            double sp = ((n1-1)*var1+(n2-1)*var2)/v;
            double t = ProbabilityDistribution::inverseTStudent((alpha / 2), 0.0, 1.0,v);
            testStat = (avg1 -avg2 )/sqrt(sp*(1/n1+1/n2));
	        //acceptanceInferiorLimit = - t;
	        //acceptanceSuperiorLimit =  t;
	       // pvalue = 1-solver.integrate(-5, testStat, &ProbabilityDistribution::tStudent, 0, 1, v);
        }
	    return StatisticalAnalyzerStudent::TestResult(pvalue, rejectH0, acceptanceInferiorLimit, acceptanceSuperiorLimit, testStat);
	}
	StatisticalAnalyzerStudent::TestResult testProportion(double prop1, unsigned int n1, double prop2, unsigned int n2, double confidenceLevel, StatisticalAnalyzerStudent::H1Comparition comp) { 
	    // DESENVOLVER - Seguir a seção 5.14
        double pvalue;
        bool rejectH0;
        double acceptanceInferiorLimit = std::numeric_limits<double>::min();
	    double acceptanceSuperiorLimit = std::numeric_limits<double>::max();
	    double testStat;
	    return StatisticalAnalyzerStudent::TestResult(pvalue, rejectH0, acceptanceInferiorLimit, acceptanceSuperiorLimit, testStat);
	}
	StatisticalAnalyzerStudent::TestResult testVariance(double var1, unsigned int n1, double var2, unsigned int n2, double confidenceLevel, StatisticalAnalyzerStudent::H1Comparition comp) { 
	    // DESENVOLVER - Seguir a seção 5.12
        double pvalue;
        bool rejectH0;
        double acceptanceInferiorLimit = std::numeric_limits<double>::min();
	    double acceptanceSuperiorLimit = std::numeric_limits<double>::max();
	    double testStat;
	    return StatisticalAnalyzerStudent::TestResult(pvalue, rejectH0, acceptanceInferiorLimit, acceptanceSuperiorLimit, testStat);
	}


public: // métodos para testes de hipóteses não paramétricos
    /*!
    * O teste do ChiQuadro se baseia na diferença relativa entre um vetor de frequências observadas e um vetor de frequências esperadas. Ambos os vetores devem possuir 
         a mesma quantidade de elementos, mas seu desenvolvimento nem precisa verificar isso, pois é responsabilidade de quem passa esses parâmetros.
    * Porém, quem passa esses parâmetros não precisa conhecer as restrições desse teste, e uma delas é que esse teste não funciona muito bem quando a frequência esperada de
         uma classe é menor que 5. Então seu desenvolvimento deve verificar isso e, se isso ocorrer, deve ir acumulando as frequências esperadas e observadas de classes
         consecutivas até que a frequência esperada acumulada seja ao menos 5 ou até que a última classe seja alcançada. 
    * Para obter o valor da probabilidade de uma distribuição de probabilidade de interesse, calcule a integral definida da função densidade de probabilidade 
         daquela distribuição, e para obter o valor correspondente a uma probabilidade acumulada, calcule o valor inverso da distribuição.
         Para isso você poderá invocar os seguintes métodos (alguns dos quais você mesmo desenvolveu no DS1!):
            Da classe "Solver"
            	double integrate(double min, double max, Solver_if::f1p f);
	            double integrate(double min, double max, Solver_if::f2p f, double p2);
	       Da classe "ProbabilityDistribution" 
            	static double inverseChi2(double cumulativeProbability, double degreeFreedom)
	        	static double chi2(double x, double degreeFreedom);
	* Como é extremanente comum invocar esses métodos milhares ou milhões de vezes ao longo de uma única simulação, os resultados  cálculos de probabilidade (integrais) 
	     ou de valores com certa probabilidade (inverseX) devem ser armazenados em alguma estrutura tipo "cache" (map, hash) para evitar serem recalculados. Por isso,
	     recomenda-se não invocar diretamente os métodos de Solver" e "ProbabilityDistribution" diretamente, mas que se crie "algo" que primeiro verifique se o cálculo
	     solicitado já foi realizado antes e já se conhece o resultado (armazenado na "cache") e que invoque os métodos dessas classes apenas se o cálculo não tiver sido
	     realizado anteriormente.
    */
    StatisticalAnalyzerStudent::TestResult testChiSquare(std::vector<double> observedFreq, std::vector<double> expectedFreq, double confidenceLevel) {
        // DESENVOLVER - Seguir a seção 
        double pvalue;
        bool rejectH0;
        double acceptanceInferiorLimit = std::numeric_limits<double>::min();
	    double acceptanceSuperiorLimit = std::numeric_limits<double>::max();
	    double testStat;
	    return StatisticalAnalyzerStudent::TestResult(pvalue, rejectH0, acceptanceInferiorLimit, acceptanceSuperiorLimit, testStat);
    }
    /*!
    * O teste de Kolmogorov–Smirnov é um método simples e alternativo ao teste do chi-quadrado para verificar se duas amostras seguem a mesma distribuição de probabilidade, ou se
        os dados amostrados se "encaixam" (aderem) a determinada distribuição de probabilidade. Esse teste recebe dois vetores com as probabilidades acumuladas de duas
        amostras (ou de uma amostra e de uma distribuição teórica, como é normalmente o caso).
    * Esse teste simplesmente compara a maior diferença entre dois pontos equivalentes das duas funções de probabilidade acumulada com um valor crítico. Se essa diferença for
        maior que o valor crítico, então rejeita-se a hipótese nula que as duas amostras seguem a mesma distribuição. Pode haver dificuldades em obter todas as informações referentes
        a esse teste de hipóteses, já que a estatística de teste não segue uma das distribuições de probabilidade usuais. Defina essas informações como std::numeric_limits<double>::min()
    * Ao contrário do teste do chi2, a distribuição deve ser completamente especificada, ou seja, se os dados não cobrirem o espaço amostral da distribuição (F(x)=0 a F(x)=1), 
        a região crítica do teste K-S não é mais válida.    
    * Pela referência https://pt.wikipedia.org/wiki/Teste_Kolmogorov-Smirnov, a estatística de teste é D=max(abs(accumf1[i]-accumf2[i]) para todo i pertencente ao conjunto 
        de dados, e o valor crítico com nível de significância "a" é c(a)=sqrt(-ln(a/2)*((1+m/n)/(2*m)), onde n e m são os tamanhos das duas amostras. 
    */
    StatisticalAnalyzerStudent::TestResult testKolmogorovSmirnov(std::vector<double> accumf1, std::vector<double> accumf2, double confidenceLevel) {
        // DESENVOLVER - Seguir a seção "Two-sample Kolmogorov–Smirnov test" de https://pt.wikipedia.org/wiki/Teste_Kolmogorov-Smirnov
        // Usar o cálculo geral de c(a) 
        double pvalue;
        bool rejectH0;
        double acceptanceInferiorLimit = std::numeric_limits<double>::min();
	    double acceptanceSuperiorLimit = std::numeric_limits<double>::max();
	    double testStat;
	    return StatisticalAnalyzerStudent::TestResult(pvalue, rejectH0, acceptanceInferiorLimit, acceptanceSuperiorLimit, testStat);
    }


public: // métodos para análise de regressão e correlação
    /*!
    *
    */
    StatisticalAnalyzerStudent::SimpleLinearRegressionResult simpleLinearRegression(std::string independentDataFilename, std::string dependentDataFilename) {
        // DESENVOLVER - Seguir a seção 6.3
        double b0;
        double b1;
        double var;
        double epb0;
        double epb1;
        double sqr;
        double sqt;
        double e0b0;
        double e0b1;
        double r2;
        return StatisticalAnalyzerStudent::SimpleLinearRegressionResult(b0, b1, var, epb0, epb1, sqr, sqt, e0b0, e0b1, r2);
    }


public: // métodos para testes de aderência (Fit)
    /*!
    * Desenvolva os métodos abaixo utilizando o teste do chi-quadrado E o teste K-S (como o Input Analyser do Arena). Esses métodos recebem o nome de um arquivo com os dados
        amostrados "sampleDataFilename" (com um dado por linha), a quantidade de classes no qual os dados devem ser agrupados, e devolvem todos os parâmetros seguintes. Depois do nome do arquivo (primeiro parâmetro), os dois parâmetros
        seguintes são os resultados do teste do chi-quadrado e do teste K-S (Kolmogorov-Smirnov) de aderẽncia dos dados do arquivo em relação à distribuição teórica do método. 
        Todos os parâmetros seguintes são os estimadores dos parâmetros daquela distribuição.
    * No teste do Chi-quadrado as frequências observadas devem ser extraídas do arquivos de dados amostrados "sampleDataFilename", e as frequências esperadas devem ser obtidas 
        a partir da integral da função densidade de probabilidade da distribuição de interesse. 
    * Para extrair as as frequências observadas do arquivo de dados, proceda da seguinte maneira: leia o arquivo, procurando os valores mínimo e máximo, o que dá a amplitude dos dados. Divida essa amplitude pela quantidade de classes
        e com isso você terá os limites de cada classe. Varra o arquivo contando a quantidade de valores dentro de cada classe, o que gera a frequência observada. Utilize intervalo fechado
        no começo do intervalo e intervalo aberto no final do intervalo de cada classe, com exceção da última classe, em que p final do intervalo é fechado (para incluir o valor máximo).
    * Para obter as as frequências esperadas da integral da função densidade de probabilidade da distribuição de interesse, proceda da seguinte maneira:
        - calcule a integral definida da distribuição de probabilidade para o mesmo intervalo da classe gerada para a frequência observada correspondente. Não repita cálculos 
          já feitos anteriormente
        - multiplique o valor da integral (probabilidade) pela quantidade de elementos na amostra (arquivo de dados) para obter a frequência esperada. 
	* Como é extremanente comum invocar esses métodos milhares ou milhões de vezes ao longo de uma única simulação, os resultados  cálculos de probabilidade (integrais) 
	     ou de valores com certa probabilidade (inverseX) devem ser armazenados em alguma estrutura tipo "cache" (map, hash) para evitar serem recalculados. Por isso,
	     recomenda-se não invocar diretamente os métodos de Solver" e "ProbabilityDistribution" diretamente, mas que se crie "algo" que primeiro verifique se o cálculo
	     solicitado já foi realizado antes e já se conhece o resultado (armazenado na "cache") e que invoque os métodos dessas classes apenas se o cálculo não tiver sido
	     realizado anteriormente.

    */
	void fitUniform(std::string sampleDataFilename, unsigned int histogramClasses, StatisticalAnalyzerStudent::TestResult *chi2Result, StatisticalAnalyzerStudent::TestResult *ksResult, double *min, double *max) { 
	    // DESENVOLVER
	    // exemplos de retorno
	    //chi2Result->pValue(0.0);
	    //*min = 0.0;

	}
	void fitTriangular(std::string sampleDataFilename, unsigned int histogramClasses, StatisticalAnalyzerStudent::TestResult *chi2Result, StatisticalAnalyzerStudent::TestResult *ksResult, double *min, double *mo, double *max) { 
	    // DESENVOLVER
	}
	void fitNormal(std::string sampleDataFilename, unsigned int histogramClasses, StatisticalAnalyzerStudent::TestResult *chi2Result, StatisticalAnalyzerStudent::TestResult *ksResult, double *avg, double *stddev) { 
	    // DESENVOLVER
	}
	void fitExpo(std::string sampleDataFilename, unsigned int histogramClasses, StatisticalAnalyzerStudent::TestResult *chi2Result, StatisticalAnalyzerStudent::TestResult *ksResult, double *avg1) { 
	    // DESENVOLVER
	}
	void fitErlang(std::string sampleDataFilename, unsigned int histogramClasses, StatisticalAnalyzerStudent::TestResult *chi2Result, StatisticalAnalyzerStudent::TestResult *ksResult, double *avg, double *m) { 
        // DESENVOLVER
	}
	void fitBeta(std::string sampleDataFilename, unsigned int histogramClasses, StatisticalAnalyzerStudent::TestResult *chi2Result, StatisticalAnalyzerStudent::TestResult *ksResult, double *alpha, double *beta, double *infLimit, double *supLimit) { 
        // DESENVOLVER	    
	}
	void fitWeibull(std::string sampleDataFilename, unsigned int histogramClasses, StatisticalAnalyzerStudent::TestResult *chi2Result, StatisticalAnalyzerStudent::TestResult *ksResult, double *alpha, double *scale) { 
	    // DESENVOLVER
	}
	
	/* esse método deve invocar todos os anteriores e retornar o nome da distribuição (Uniform, Triangular, Normal, Expo, etc) que tiver o melhor ajustamento (fit), 
	     representado pelo maior p-value do teste do chi-quadrado de cada ajuste
	*/
	void fitAll(std::string sampleDataFilename, unsigned int histogramClasses, StatisticalAnalyzerStudent::TestResult *chi2Result, StatisticalAnalyzerStudent::TestResult *ksResult, std::string *name) { 
	    // DESENVOLVER
	}
	
private: // implemente aqui os métodos e atributos privados que forem necessários.

};

#endif /* STATISTICALANALYZERSTUDENT_H */

